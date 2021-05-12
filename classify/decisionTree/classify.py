import sklearn
import numpy as np
import time
from sklearn import tree
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score

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

def random_forest():
	print("---------------------------")
	start_time = time.time()

	num_trees = 51

	clf = RandomForestClassifier(n_estimators=num_trees, n_jobs=-1, max_features=0.5)
	clf.fit(train_X, train_Y)
	
	preds = clf.predict(test_X)
	correct = (preds == test_Y).sum()
	accuracy = correct / len(test_Y)

	end_time = time.time()
	print("Random forest statistics with " + str(num_trees) + " trees:")
	print(str(correct) + "/" + str(len(preds)) + " correct. Acc:", accuracy)
	print(classification_report(test_Y, preds))
	print("Elapsed time:", end_time - start_time)

random_forest()

