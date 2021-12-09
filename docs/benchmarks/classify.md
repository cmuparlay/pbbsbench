---
title: Classify (CLAS)
---

# Classify (CLAS) 

This supervised-learning benchmark chooses labels for a set of
unlabeled feature vectors after being trained on a set of labeled
feature vectors.

The input consists of a training set and a test set.  The training set
consists of a sequence of n feature vectors each of length k and each
with an associated label.
Each column (position in the feature vector) is either discrete or
continuous.  An indicator vector of length k is used to specify the
type of each column.
The label is always discrete.
For both discrete and continuous vectors, the values in both the
feature vectors and labels have to be in the range [0-256).
The test set is a sequence of m feature vectors of each of length k, and without labels.

The output is a sequence of predictions for the labels for each of the m test vectors.

The data sets include "the ground truth" for each of the labels of the test set, which can then
be compared to the predicted labels to report accuracy.

We note that there is no correct answer for this benchmark, so results should be measured as
a point or set of points on an accuracy/time grid.    Accuracy is simply the percent correct.

The testing harness reports the percent correct.

### Default Input Distributions

The test distributions are two datasets take from the [UC Irvine Machine Learning repository](https://archive.ics.uci.edu/ml/index.php).  The data has been normalized so that the values are in the specified range.  This was done using scalar quantization so that approximately equal number of items fall in each bucket.   It is not necessarily a linear transform.   The instances are:

- [`covtype.data`](https://archive.ics.uci.edu/ml/datasets/covertype) : this is a data set used to predict the cover type of
  forests from cartographic variables (e.g. altitude, soil types,
  slope directions, etc.).

- [`kddcup.data`](https://archive.ics.uci.edu/ml/datasets/kdd+cup+1999+data) : this data set is from the kdd cup competition
  in 1999. The goal is to predict whether network packets are malicious
  or not.

All instances `<instance>` include three files `<instance>.train`,
`<instance>.test`, and `<instance>.labels`.

### Input and Output File Formats 

All files are in ascii.

The `.train` file consist of a first line which is a comma separated
sequence of length k each containing either `i` (continuous) or `d`
(discrete).
The next n lines are the feature vectors and the label, each consisting of comma
separated sequence of k + 1 integers where the label is at the end.

The `.test` file is the same format except it is missing the labels.

The `.labels` file is a comma separated sequence of integers of length
n. The output should be in this format.



The input is an ascii string containing the documents.   The output is
as ascii string as described above. 
