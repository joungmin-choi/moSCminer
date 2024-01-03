#####################################
### Specify your input parameters ###
#####################################
num_omics = 3 # Number of your omics
dirname = "./example_dataset/"
resDir = "./results/"
#####################################

#from __future__ import print_function
import tensorflow as tf
import numpy as np
import pandas as pd
from pandas import DataFrame as df
import os
import sys
from tensorflow.contrib.layers.python.layers import batch_norm as batch_norm
from sklearn.model_selection import train_test_split

# SET ENV
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
config = tf.ConfigProto()
config.intra_op_parallelism_threads = 44
config.inter_op_parallelism_threads = 44
config.gpu_options.allow_growth=True

# FUNCTIONS
def fc_bn(_x, _output, _phase, _scope):
	with tf.variable_scope(_scope):
		h1 = tf.contrib.layers.fully_connected(_x, _output, activation_fn=None, scope='dense', weights_initializer=tf.contrib.layers.variance_scaling_initializer(), weights_regularizer = tf.contrib.layers.l2_regularizer(0.01))
		h2 = tf.layers.batch_normalization(h1, fused=True, center=True, scale=True, training=_phase, name='bn')
		return h2

try :
	os.mkdir(resDir)
except :
	print("Exist!")


x_omics_train_df_list = []
x_omics_test_df_list = []
for i in range(1, num_omics+1) :
	temp_train = pd.read_csv(os.path.join(dirname, "final_" + str(i) + "_train.csv"), index_col = 0)
	temp_test = pd.read_csv(os.path.join(dirname, "final_" + str(i) + "_test.csv"), index_col = 0)
	x_omics_train_df_list.append(temp_train)
	x_omics_test_df_list.append(temp_test)

y_data = pd.read_csv(dirname + "final_celltype_onehot_train.csv", index_col = 0)
y_test = pd.read_csv(dirname + "final_celltype_onehot_test.csv", index_col = 0)

num_feature_list = []

for j in range(num_omics) :
	num_feature_list.append(len(x_omics_train_df_list[j].columns))

n_sample = len(x_omics_train_df_list[0])
n_classes = len(y_data.columns)

# PLACEHOLDER

tf_raw_X_list = []
for j in range(num_omics) :
	tf_raw_X_list.append(tf.placeholder(tf.float32, [None, num_feature_list[j]]))

tf_raw_Y = tf.placeholder(tf.float32, [None, n_classes])
keep_prob = tf.placeholder(tf.float32)
phase = tf.placeholder(tf.bool, name='phase')
handle = tf.placeholder(tf.string, shape=[])
noise_r = tf.placeholder(tf.float32, shape=[])

tf_raw_X_list.append(tf_raw_Y)
tf_raw_X_tuple = tuple(tf_raw_X_list)

# PARAMETERS
if n_sample > 128 :
	batch_size = 128 
else :
	batch_size = n_sample

ep_repeat = 10 
prefetch_size = batch_size * 2
test_data_size = len(x_omics_test_df_list[0])

learn_rate_sm = 1e-3
keep_rate_sm = 0.5

ensemble_model_num = 1
train_sm_eps = 1000

# DATASET & ITERATOR
dataset_train = tf.data.Dataset.from_tensor_slices(tf_raw_X_tuple)
dataset_train = dataset_train.shuffle(buffer_size=batch_size * 2)
dataset_train = dataset_train.repeat(ep_repeat).batch(batch_size).prefetch(prefetch_size)
iterator_train = dataset_train.make_initializable_iterator()

dataset_test = tf.data.Dataset.from_tensor_slices(tf_raw_X_tuple)
dataset_test = dataset_test.batch(test_data_size)
iterator_test = dataset_test.make_initializable_iterator()

iter = tf.data.Iterator.from_string_handle(handle, dataset_train.output_types, dataset_train.output_shapes)

*tf_X_list, tf_Y = iter.get_next()

# MODEL STRUCTURE
n_sm_h2 = 200 
n_sm_out = n_classes

n_embedding = 128 #8 #16 #64
n_prj = 64 #8 #32

# MODEL FUNCTIONS
def attention(_code, n_features, _keep_prob, _phase, omic_name):
	embedding_vector = np.random.rand(n_features, n_embedding).astype(np.float32)
	f_e = _code[:,:,None] * embedding_vector
	f_x = tf.nn.dropout(tf.nn.tanh(fc_bn(f_e, n_prj, _phase, "f_x_" + omic_name)), _keep_prob)
	f_a_1 = tf.nn.dropout(tf.nn.tanh(fc_bn(f_e, n_features, _phase, "f_a_1_" + omic_name)), _keep_prob)
	f_a_2 = fc_bn(f_a_1, 1, _phase, "f_a_2_" + omic_name)
	f_a_2 = tf.reshape(f_a_2, shape = [-1,n_features])
	importance = tf.nn.softmax(f_a_2)
	new_representation = tf.reshape(importance, [-1,n_features,1]) * f_x
	new_representation = tf.reduce_sum(new_representation, 1)
	return new_representation, importance

def softmax(rep_omics_list, num_omics, _keep_prob, _phase):
	rep_concated = rep_omics_list[0]
	for j in range(1, num_omics) :
		rep_concated = tf.concat([rep_concated, rep_omics_list[j]], axis = 1)
	fc2 = tf.nn.dropout(tf.nn.elu(fc_bn(rep_concated, n_sm_h2, _phase, "fc2")), _keep_prob)
	sm = tf.nn.softmax(fc_bn(fc2, n_sm_out, _phase, "sm"))
	return sm

# MODEL
rep_omics_list = []
attn_omics_list = []
for j in range(num_omics) :
	temp_rep, temp_attn = attention(tf_X_list[j], num_feature_list[j], keep_prob, phase, "omics" + str(j))
	rep_omics_list.append(temp_rep)
	attn_omics_list.append(temp_attn)


sm_out = softmax(rep_omics_list, num_omics, keep_prob, phase)
print ("MODEL READY") 

# DEFINE LOSS AND OPTIMIZER
update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
with tf.control_dependencies(update_ops) :
	sm_cost = tf.reduce_mean(-tf.reduce_sum(tf_Y * tf.log(sm_out + 1e-10), axis = 1))
	train_op_sm = tf.train.AdamOptimizer(learning_rate=learn_rate_sm).minimize(sm_cost)


# ACCURACY
pred = tf.argmax(sm_out, 1)
label = tf.argmax(tf_Y, 1)
correct_pred = tf.equal(pred, label)
accuracy = tf.reduce_mean(tf.cast(correct_pred, tf.float32))
_accuracy = tf.Variable(0)

print ("FUNCTIONS READY") 

# START SESSION
sess = tf.Session(config=config)
handle_train = sess.run(iterator_train.string_handle())
handle_test = sess.run(iterator_test.string_handle())

print ("START OPTIMIZATION & TESTING")
for model_num in range(ensemble_model_num):
	sess.run(tf.global_variables_initializer())	

	# SET OPS & FEED_DICT
	sm_ops = [sm_cost, train_op_sm, accuracy]
	sm_ops.extend(attn_omics_list)

	sm_feed_dict = {handle: handle_train, keep_prob: keep_rate_sm, phase: True}
	for temp_ep, meta_step, temp_ops, temp_feed_dict in zip([train_sm_eps], ["_sm"], [sm_ops], [sm_feed_dict]):	
		for ep in xange(temp_ep//ep_repeat):
			feed_dict_train = {}
			for j in range(num_omics) :
				feed_dict_train[tf_raw_X_list[j]] = x_omics_train_df_list[j]
			feed_dict_train[tf_raw_Y] = y_data

			sess.run(iterator_train.initializer, feed_dict = feed_dict_train)

			# REPEAT NUMBER OF EPS WITHOUT BREAK. SET BY ep_repeat
			while True: 
				try:
					cur_cost_val, _, cur_accuracy, *cur_attn_list = sess.run(temp_ops, feed_dict = temp_feed_dict)
				except tf.errors.OutOfRangeError:
					break

			# EXECUTED PER ep_repeat
			print("Ep:%04d" % (ep*ep_repeat), ", Cost:" + str(cur_cost_val) + ", Train_batch_accr:" + str(cur_accuracy))
		
		feed_dict_test = {}
		for j in range(num_omics) :
			feed_dict_test[tf_raw_X_list[j]] = x_omics_test_df_list[j]
		feed_dict_test[tf_raw_Y] = y_test

		sess.run(iterator_test.initializer, feed_dict = feed_dict_test)
		cur_pred = sess.run([pred], feed_dict = {handle: handle_test, keep_prob: 1.0, phase: False})


celltype_list = y_data.columns.tolist()
cur_pred_str_list = []
for cell in cur_pred[0] :
	cur_pred_str_list.append(celltype_list[cell])

cur_pred_df = pd.DataFrame({'SampleName' : x_omics_test_df_list[0].index.tolist(), 'cellType' : cur_pred_str_list})
cur_pred_df.to_csv(resDir + "test_data_celltype_prediction_results.csv", mode = "w", index = False)


for i in range(len(cur_attn_list)) :
	cur_attn_omics_df = pd.DataFrame(cur_attn_list[i])
	cur_attn_omics_df.columns = x_omics_train_df_list[i].columns
	cur_attn_omics_df_rank = pd.DataFrame(cur_attn_omics_df.sum())
	cur_attn_omics_df_rank.sort_values(by = 0, ascending = False, inplace = True)
	cur_attn_omics_df_rank.reset_index(inplace = True, drop = False)
	cur_attn_omics_df_rank['rank'] = list(range(1,len(cur_attn_omics_df_rank)+1))
	del cur_attn_omics_df_rank[0]
	cur_attn_omics_df_rank.to_csv(resDir +"omics" + str(i+1) + "_importance_rank.csv", mode = "w", index = False)
