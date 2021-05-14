import sklearn
import numpy as np
import time
from sklearn import tree
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import classification_report, accuracy_score
from sklearn.neural_network import MLPClassifier

train_file_path = "../sequenceData/data/covtype.data.train"
test_file_path = "../sequenceData/data/covtype.data.test"

def load_data(file_path):
	data = np.genfromtxt(file_path, delimiter=",")
	data = data[1:] # Remove i, d labels
	X = data[:, :(data.shape[1] - 1)] #last column is labels
	Y = data[:, data.shape[1] - 1]
	return X, Y

train_X, train_Y = load_data(train_file_path)
test_X, test_Y = load_data(test_file_path)

print("Loaded data")

def decision_tree():
	print("-------------------")
	start_time = time.time()

	clf = tree.DecisionTreeClassifier()
	clf = clf.fit(train_X, train_Y)

	preds = clf.predict(test_X)
	correct = (preds == test_Y).sum()
	accuracy = correct / len(test_Y)

	end_time = time.time()
	print("Decision tree statistics:")
	print(str(correct) + "/" + str(len(preds)) + " correct. Acc:", accuracy)
	print(classification_report(test_Y, preds))
	print("Elapsed time:", end_time - start_time)

decision_tree()

def random_forest(num_trees):
	print("---------------------------")
	start_time = time.time()

	clf = RandomForestClassifier(n_estimators=num_trees, n_jobs=-1, max_features=0.7)
	clf.fit(train_X, train_Y)
	
	preds = clf.predict(test_X)
	correct = (preds == test_Y).sum()
	accuracy = correct / len(test_Y)

	end_time = time.time()
	print("Random forest statistics with " + str(num_trees) + " trees:")
	print(str(correct) + "/" + str(len(preds)) + " correct. Acc:", accuracy)
	print(classification_report(test_Y, preds))
	print("Elapsed time:", end_time - start_time)

random_forest(81)
random_forest(243)

#Note: very slow and occasionally crashes due to having many jobs
def knn(k):
	print("---------------------------")
	start_time = time.time()

	clf = KNeighborsClassifier(n_neighbors=k, n_jobs=72)
	clf.fit(train_X, train_Y)
	
	preds = clf.predict(test_X)
	correct = (preds == test_Y).sum()
	accuracy = correct / len(test_Y)

	end_time = time.time()
	print("k nearest neighbors with k=" + str(k) + ":")
	print(str(correct) + "/" + str(len(preds)) + " correct. Acc:", accuracy)
	print(classification_report(test_Y, preds))
	print("Elapsed time:", end_time - start_time)

def extra_forest(num_trees, max_f, to_bootstrap):
	print("---------------------------")
	start_time = time.time()

	clf = ExtraTreesClassifier(n_estimators=num_trees, n_jobs=-1, bootstrap=to_bootstrap, max_features=max_f)
	clf.fit(train_X, train_Y)
	
	preds = clf.predict(test_X)
	correct = (preds == test_Y).sum()
	accuracy = correct / len(test_Y)

	end_time = time.time()
	print("extra forest with k=" + str(num_trees))
	print(str(correct) + "/" + str(len(preds)) + " correct. Acc:", accuracy)
	print(classification_report(test_Y, preds))
	print("Elapsed time:", end_time - start_time)

extra_forest(27, None, False)
extra_forest(81, None, False)
extra_forest(243, None, False)
#extra_forest(27, 0.8, False)
#extra_forest(81, 0.6, False)
#extra_forest(81, 0.4, False)
#extra_forest(81, 0.2, False)
#extra_forest(243, "sqrt", False)
#extra_forest(243, "log2", False)
#extra_forest(27, None, True)
#extra_forest(81, None, True)
#extra_forest(243, None, True)

def mlp(hidden_sizes=(200, 200), lr=0.001, ):
	print("---------------------------")
	start_time = time.time()

	clf = MLPClassifier(hidden_layer_sizes=hidden_sizes, learning_rate_init=lr, max_iter=50, verbose=True)
	clf.fit(train_X, train_Y)
	
	preds = clf.predict(test_X)
	correct = (preds == test_Y).sum()
	accuracy = correct / len(test_Y)

	end_time = time.time()
	print("NN with sizes=" + str(hidden_sizes))
	print(str(correct) + "/" + str(len(preds)) + " correct. Acc:", accuracy)
	print(classification_report(test_Y, preds))
	print("Elapsed time:", end_time - start_time)

#mlp()
