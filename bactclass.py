"""!
BactClass: Machine Learning Classifier for Bactome

Date created: 25th July 2020

License: GNU General Public License version 3 for academic or 
not-for-profit use only


Bactome package is free software: you can redistribute it and/or 
modify it under the terms of the GNU General Public License as 
published by the Free Software Foundation, either version 3 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
import pickle

import pandas as pd
from sklearn.model_selection import train_test_split

try: 
    import fire
except ImportError:
    subprocess.check_call([sys.executable, '-m', 'pip', 
                           'install', 'fire'])
    import fire

try: 
    import joblib
except ImportError:
    subprocess.check_call([sys.executable, '-m', 'pip', 
                           'install', 'joblib'])
    import joblib

def readData(filename, training=True, label="Class"):
    """!
    Internal function - Read the data file into a Pandas dataframe.

    @param filename String: Path to CSV data file.
    @param training Boolean: Flag to indicate whether the file consists of data for training. Default = True
    @param label String: Column (field) name in the data file to indicate the class label. Default = Class
    """
    data = pd.read_csv(filename)
    if training:
        X = data.drop(label, axis=1)
        Y = data[label]
        return X, Y
    else: 
        return data
    
def saveModel(filename, filetype, classifier):
    """!
    Internal function - Write out the classifier object into a pickle file.

    @param filename String: Path to write out the generated classifier.
    @param filetype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib".
    @param classifier Object: Classifier object
    """
    if filetype.lower() == "pickle":
        f = open(filename, "wb")
        pickle.dump(classifier, f, pickle.HIGHEST_PROTOCOL)
        f.close()
    elif filetype.lower() == "joblib":
        joblib.dump(classifier, filename)

def loadModel(filename, filetype):
    """!
    Internal function - Load the classifier from file.

    @param filename String: Path to the generated classifier to be loaded.
    @param filetype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib". 
    """
    if filetype.lower() == "pickle":
        f = open(filename, "rb")
        classifier = pickle.load(f)
        f.close()
        return classifier
    elif filetype.lower() == "joblib":
        classifier = joblib.load(filename)
        return classifier

def showClassifierParameters(classifier):
    """!
    Internal function - Print out parameters of the classifier.

    @param classifier Object: Classifier object
    """
    param = classifier.get_params()
    print("------------- Classifier Parameters -----------------")
    for k in param:
        print("%s = %s" % (k, param[k]))
    print("------------- End of Classifier Parameters ----------")
    print("")
    
def showConfusionMatrix(Y_test, Y_pred):
    """!
    Internal function - Print out confusion matrix.

    @param Y_test Array: Actual classification from training data in NumPy Array.
    @param Y_pred Array: Predicted classification from classifier in NumPy Array.
    """
    from sklearn.metrics import confusion_matrix
    matrix = confusion_matrix(Y_test, Y_pred)
    print("------------- Confusion Matrix ----------------------")
    for row in matrix:
        print(', '.join([str(x) for x in row]))
    print("------------- End of Confusion Matrix ---------------")
    print("")
    
def showClassificationReport(Y_test, Y_pred):
    """!
    Internal function - Print out classification report (precision, recall, F1-score, and Accuracy).

    @param Y_test Array: Actual classification from training data in NumPy Array.
    @param Y_pred Array: Predicted classification from classifier in NumPy Array.
    """
    from sklearn.metrics import classification_report
    print("------------- Classification Report -----------------")
    print(classification_report(Y_test, Y_pred))
    print("------------- End of Classification Report ----------")
    print("")

def cross_validate(classifier, X, Y, fold):
    """!
    Internal function - Cross validation of classifier

    @param classifier Object: Generated classifier object.
    @param X Array: Array of independent variables.
    @param Y Array: Array of dependent variable.
    @param fold Integer: Number of folds for cross validation
    """
    from sklearn import metrics
    from sklearn.model_selection import cross_val_score
    fold = int(fold)
    print("------------ Cross Validation Report ----------------")
    print("%s fold cross validation" % str(fold))
    print("")
    scores = cross_val_score(classifier, X, Y, cv=fold, scoring="precision")
    print("Precision:                     %0.3f (sigma = %0.4f)" % (scores.mean(), scores.std()))
    scores = cross_val_score(classifier, X, Y, cv=fold, scoring="recall")
    print("Recall:                        %0.3f (sigma = %0.4f)" % (scores.mean(), scores.std()))
    scores = cross_val_score(classifier, X, Y, cv=fold, scoring="f1")
    print("F1 Score:                      %0.3f (sigma = %0.4f)" % (scores.mean(), scores.std()))
    scores = cross_val_score(classifier, X, Y, cv=fold, scoring="accuracy")
    print("Accuracy:                      %0.3f (sigma = %0.4f)" % (scores.mean(), scores.std()))
    scores = cross_val_score(classifier, X, Y, cv=fold, scoring="roc_auc")
    print("Area Under ROC Curve (AUC):    %0.3f (sigma = %0.4f)" % (scores.mean(), scores.std()))
    scores = cross_val_score(classifier, X, Y, cv=fold, scoring="jaccard")
    print("Jaccard Similarity:            %0.3f (sigma = %0.4f)" % (scores.mean(), scores.std()))
    scores = cross_val_score(classifier, X, Y, cv=fold, scoring="normalized_mutual_info_score")
    print("Normalized Mutual Information: %0.3f (sigma = %0.4f)" % (scores.mean(), scores.std()))
    scores = cross_val_score(classifier, X, Y, cv=fold, scoring="v_measure_score")
    print("V-measure:                     %0.3f (sigma = %0.4f)" % (scores.mean(), scores.std()))
    scores = cross_val_score(classifier, X, Y, cv=fold, scoring="neg_mean_squared_error")
    print("Mean Square Error (MSE):       %0.3f (sigma = %0.4f)" % (-1*scores.mean(), scores.std()))
    scores = cross_val_score(classifier, X, Y, cv=fold, scoring="r2")
    print("R-square:                      %0.3f (sigma = %0.4f)" % (scores.mean(), scores.std()))
    print("--------- End of Cross Validation Report ------------")

def process_classifier(datafile, label, classifier, oclass, otype, 
                        classparam, confusion, classreport, cross_validation):
    """!
    Internal function - Process (train, test, and save) classifier.

    @param datafile String: Path to CSV data file used to generate SVM.
    @param label String: Column (field) name in the data file to indicate the class label.
    @param oclass String: Path to write out the generated classifier. 
    @param otype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib".
    @param classifier Object: Generated classifier object.
    @param classparam Boolean: Flag to indicate whether to print out SVM parameters.
    @param confusion Boolean: Flag to indicate whether to print out confusion matrix. 
    @param classreport Boolean: Flag to indicate whether to print out classification report.
    @param cross_validation Integer: Number of cross validation (if any) to perform. If less than 2, cross validation will not be carried out.
    """
    X, Y = readData(datafile, True, label)
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size = 0.20)
    classifier.fit(X_train, Y_train)
    Y_pred = classifier.predict(X_test)
    saveModel(oclass, otype, classifier)
    if classparam:
        showClassifierParameters(classifier)
    if confusion:
        showConfusionMatrix(Y_test, Y_pred)
    if classreport:
        showClassificationReport(Y_test, Y_pred)
    if cross_validation > 1:
        cross_validate(classifier, X, Y, int(cross_validation))

def useScikitClassifier(classifier_type, datafile, classfile, classtype, resultfile):
    """!
    Internal function - Use a previously generated SciKit-Learn classifier  
    to classify data.

    @param classifier_type String: Type of classifier.
    @param datafile String: Path to CSV file containing data to be classified.
    @param classfile String: Path to the generated classifier.
    @param classtype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib".
    @param resultfile String: Path to write out the classified results.
    """
    taskText = {"ANN": "Classifying using Artificial Neural Network (ANN)",
                "DT": "Classifying using Decision Tree (DT)",
                "SVM": "Classifying using Support Vector Machine (SVM)",
                "GNB":"Classifying using Gaussian Naive Bayes (GNB)",
                "BNB":"Classifying using Bernoulli Naive Bayes (BNB)",
                "MNB":"Classifying using Multinomial Naive Bayes (MNB)",
                "CNB":"Classifying using Complementary Naive Bayes (CNB)"}
    print("")
    print("Task: %s" % taskText[classifier_type])
    print("Parameters:")
    print("    Data File = " + str(datafile))
    print("    Classifier File = " + str(classfile))
    print("    Classifier File Type = " + str(classtype))
    print("    Result File = " + str(resultfile))
    print("")
    classifier = loadModel(classfile, classtype)
    data = readData(datafile, False)
    data["Prediction"] = classifier.predict(data)
    data.to_csv(resultfile, index=False)

def recycle(infile, intype, outfile, outtype):
    """!
    Function to read in and write out a classifier. This can be used to update the serialization protocol or to change the change the type of serialization; such as, from pickle to joblib.

    Usage:
        
        python bactclass.py recycle --infile=classifier_ANN.pickle --intype=pickle --outfile=classifier_ANN.joblib --outtype=joblib

    @param infile String: Path to the input classifier.
    @param intype String: Type of file of the input classifier. Allowable types are "pickle" and "joblib".
    @param outfile String: Path to the output classifier.
    @param outtype String: Type of file of the output classifier. Allowable types are "pickle" and "joblib".
    """
    print("")
    print("Task: Read and Write (Recycle) a Classifier")
    print("Parameters:")
    print("    Input Classifier File = " + str(infile))
    print("    Input Classifier File Type = " + str(intype))
    print("    Output Classifier File = " + str(outfile))
    print("    Output Classifier File Type = " + str(outtype))
    print("")
    classifier = loadModel(infile, intype)
    saveModel(outfile, outtype, classifier)

def generateANN(datafile, label, 
                oclass="classifier_ANN.pickle", 
                otype="pickle",
                hidden_layer_sizes="100",
                activation="relu",
                solver="adam",
                learning_rate="constant",
                learning_rate_init=0.001,
                power_t=0.5,
                max_iteration=200,
                shuffle=True,
                tolerance=0.001,
                momentum=0.9,
                nesterovs_momentum=True,
                beta_1=0.9,
                beta_2=0.999,
                epsilon=1e-8,
                n_iter_no_change=10,
                verbose=False,
                classparam=True, 
                confusion=True, 
                classreport=True,
                cross_validation=5):
    """!
    Function to generate an artificial neural network (multi-layer perceptron classifier) from given data. 

    Usage:
        
        python bactclass.py genANN --datafile=classifier_train.csv --label=Class --oclass=classifier_ANN.pickle --otype=pickle --hidden_layer_sizes=100 --activation=relu --solver=adam --learning_rate=constant --learning_rate_init=0.001 --power_t=0.5 --max_iteration=200 --shuffle=True --tolerance=0.001 --momentum=0.9 --nesterovs_momentum=True --beta_1=0.9 --beta_2=0.999 --epsilon=1e-8 --n_iter_no_change=10 --verbose=False --classparam=True --confusion=True --classreport=True --cross_validation=5

    @param datafile String: Path to CSV data file used to generate SVM.
    @param label String: Column (field) name in the data file to indicate the class label.
    @param oclass String: Path to write out the generated classifier. Default = classifier_ANN.pickle
    @param otype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib". Default = pickle
    @param hidden_layer_sizes String: Defines the number of hidden layers and the number of nodes per hidden layer in semicolon-delimited format. For example, "100" represents 1 hidden layer of 100 nodes; whereas, "100;30" represents 2 hidden layers where the first hidden layer contains 100 nodes while the second hidden layer contains 30 nodes. Default = 100
    @param activation String: Activation function of the hidden layer. Allowable types are "identity", f(x) = x; "logistic" (logistic sigmoid function), f(x) = 1 / (1 + exp(-x)); "tanh" (hyperbolic tangeant function), f(x) = tanh(x); and "relu" (rectified linear unit function),f(x) = max(0, x). Default = relu
    @param solver String: Solver for weight optimization. Allowable types are "lbfgs" (is an )optimizer in the family of quasi-Newton methods), "sgd" (stochastic gradient descent), and "adam" (stochastic gradient-based optimizer proposed by Kingma and Ba, 2015. 3rd International Conference on Learning Representations). Default = adam
    @param learning_rate String: Learning rate schedule for weight updates, which is used only when solver is "sgd". Allowable types are "constant" (constant learning rate given by "learning_rate_init"), "invscaling" (gradually decreases the learning rate at each time step "t" using an inverse scaling exponent of "power_t" where ffective_learning_rate = learning_rate_init / pow(t, power_t)), and "adaptive" (keeps the learning rate constant to "learning_rate_init" as long as training loss keeps decreasing). Default = constant
    @param learning_rate_init Float: Initial learning rate used to control the step-size in updating the weights. This is only used when solver is "sgd or "adam". Default = 0.001
    @param power_t Float: Exponent for inverse scaling learning rate, which is used when solver is "sgd" for updating effective learning rate when the learning_rate is set to "invscaling". Default = 0.5
    @param max_iteration Integer: Maximum number of iterations where the solver iterates until convergence (determined by "tolerance") or reaching maximum iterations. For stochastic solvers ("sgd" or "adam"), note that this determines the number of epochs (how many times each data point will be used) rather than the number of gradient steps.
    @param shuffle Boolean: Determines whether to shuffle samples in each iteration. Only used when solver is "sgd" or "adam". Default = True
    @param tolerance Float: Tolerance for the optimization. When the loss or score is not improving by at least tolerance for n_iter_no_change consecutive iterations, unless learning_rate is set to "adaptive", convergence is considered to be reached and training stops. Default = 0.001
    @param momentum Float: Momentum, between 0 and 1, for gradient descent update (solver is "sgd"). Default = 0.9
    @param nesterovs_momentum Boolean: Flag to whether to use Nesterov's momentum when solver is "sgd" and momentum > 0. Default = True
    @param beta_1 float: Exponential decay rate for estimates of first moment vector, between 0 (inclusive) and 1 (exclusive), in "adam" solver. Default = 0.9
    @param beta_2 float: Exponential decay rate for estimates of second moment vector, between 0 (inclusive) and 1 (exclusive), in "adam" solver. Default = 0.999
    @param epsilon float: Value for numerical stability in "adam". Default = 1e-8
    @param n_iter_no_change Integer: Maximum number of epochs to not meet tolerance improvement when solver is "sgd" or "adam". Default = 10
    @param verbose Boolean: Flag to indicate whether to print progress messages. Default = False
    @param classparam Boolean: Flag to indicate whether to print out SVM parameters. Default = True
    @param confusion Boolean: Flag to indicate whether to print out confusion matrix. Default = True
    @param classreport Boolean: Flag to indicate whether to print out classification report. Default = True
    @param cross_validation Integer: Number of cross validation (if any) to perform. If less than 2, cross validation will not be carried out. Default = 5
    """
    from sklearn.neural_network import MLPClassifier
    print("")
    print("Task: Generate Artificial Neural Network (ANN) Multi-Layer Perceptron Classifier")
    print("Parameters:")
    print("    Data File = " + str(datafile))
    print("    Classification Label = " + str(label))
    print("    Hidden Layer Sizes = " + str(hidden_layer_sizes))
    print("    Activation Function = " + str(activation))
    print("    Solver Function = " + str(solver))
    print("    Type of Learning Rate = " + str(learning_rate))
    print("    Initial Learning Rate (learning_rate_init) = " + str(learning_rate_init))
    print("    Inverse Scaling Learning Rate Exponent (power_t) = " + str(power_t))
    print("    Momentum = " + str(momentum))
    print("    Nesterov's Momentum = " + str(nesterovs_momentum))
    print("    Exponential decay rate of first moment vector (beta_1) = " + str(beta_1))
    print("    Exponential decay rate of second moment vector (beta_2) = " + str(beta_2))
    print("    Numerical Stability (epsilon) = " + str(beta_2))
    print("    Maximum Iteration = " + str(int(max_iteration)))
    print("    Maximum Iteration with No Change (n_iter_no_change) = " + str(n_iter_no_change))
    print("    Tolerance = " + str(tolerance))
    print("    Shuffle = " + str(shuffle))
    print("    Verbose = " + str(verbose))
    print("    Classifier File = " + str(oclass))
    print("    Classifier File Type = " + str(otype))
    print("")
    hidden_layer_sizes = str(hidden_layer_sizes)
    hidden_layer_sizes = hidden_layer_sizes.split(";")
    hidden_layer_sizes = [int(x.strip()) for x in hidden_layer_sizes]
    classifier = MLPClassifier(hidden_layer_sizes=tuple(hidden_layer_sizes),
                               activation=str(activation),
                               solver=str(solver),
                               learning_rate=str(learning_rate),
                               learning_rate_init=float(learning_rate_init),
                               power_t=float(power_t),
                               random_state=1, 
                               max_iter=int(max_iteration),
                               shuffle=shuffle,
                               tol=float(tolerance),
                               verbose=verbose,
                               momentum=float(momentum),
                               nesterovs_momentum=nesterovs_momentum,
                               beta_1=float(beta_1),
                               beta_2=float(beta_2),
                               epsilon=float(epsilon),
                               n_iter_no_change=int(n_iter_no_change))
    process_classifier(datafile, label, classifier, oclass, otype, 
                       classparam, confusion, classreport, cross_validation)
    print("===================== ANN Generated =====================")

def useANN(datafile, classfile, classtype, resultfile):
    """!
    Function to use a previously generated artificial neural network (ANN) 
    to classify data.

    Usage:

        python bactclass.py useANN --datafile=classifier_use.csv --classfile=classifier_ANN.pickle --classtype=pickle --resultfile=classifier_result.csv
    
    @param datafile String: Path to CSV file containing data to be classified.
    @param classfile String: Path to the generated classifier.
    @param classtype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib".
    @param resultfile String: Path to write out the classified results.
    """
    useScikitClassifier("ANN", datafile, classfile, classtype, resultfile)

def generateSVM(datafile, label, 
                oclass="classifier_SVM.pickle", 
                otype="pickle",
                kernel="linear", 
                degree=3,
                gamma="scale",
                coef0=0.0,
                decision_function_shape="ovr",
                tolerance=0.001,
                max_iteration=-1,
                break_ties=False,
                classparam=True, 
                confusion=True, 
                classreport=True,
                cross_validation=5):
    """!
    Function to generate a support vector machine from given data. 

    Usage:
        
        python bactclass.py genSVM --datafile=classifier_train.csv --label=Class --oclass=classifier_SVM.pickle --otype=pickle --kernel=linear --degree=3 --gamma=scale --coef0=0.0 --decision_function_shape=ovr --tolerance=0.001 --max_iteration=-1 --break_ties=False --classparam=True --confusion=True --classreport=True --cross_validation=5

    @param datafile String: Path to CSV data file used to generate SVM.
    @param label String: Column (field) name in the data file to indicate the class label.
    @param oclass String: Path to write out the generated classifier. Default = classifier_SVM.pickle
    @param otype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib". Default = pickle
    @param kernel String: Type of SVM kernel. Acceptable values are linear, poly, rbf, sigmoid. Default = linear
    @param degree Integer: Polynomial degree, only used in polynomial (poly) kernel. Default = 3
    @param gamma String: Coefficient for "rbf", "poly" and "sigmoid" kernels. Allowable types are "scale" and "float". Default = scale
    @param coef0 float: Independent term used in "poly" and "sigmoid" kernels. Default = 0.0
    @param decision_function_shape String: Uses one-vs-rest (ovr) decision function of shape or one-vs-one ("ovo") decision function. However, one-vs-one ("ovo") is always used as multi-class strategy. The parameter is ignored for binary classification. Default = ovr
    @param tolerance float: Tolerance for stopping criterion. Default = 0.001
    @param max_iteration Integer: Hard limit on iterations within solver, or -1 for no limit. Default = -1
    @param break_ties Boolean: If True, decision_function_shape is "ovr, and number of classes > 2, prediction function will break ties according to the confidence values of decision_function; otherwise the first class among the tied classes is returned. Default = False
    @param classparam Boolean: Flag to indicate whether to print out SVM parameters. Default = True
    @param confusion Boolean: Flag to indicate whether to print out confusion matrix. Default = True
    @param classreport Boolean: Flag to indicate whether to print out classification report. Default = True
    @param cross_validation Integer: Number of cross validation (if any) to perform. If less than 2, cross validation will not be carried out. Default = 5
    """
    from sklearn.svm import SVC
    print("")
    print("Task: Generate Support Vector Machine (SVM) Classifier")
    print("Parameters:")
    print("    Data File = " + str(datafile))
    print("    Classification Label = " + str(label))
    print("    SVM Kernel Type = " + str(kernel))
    print("    Polynomial Degree = " + str(int(degree)))
    print("    Kernel Coefficient (gamma) = " + str(gamma))
    print("    Independent Coefficient (coef0) = " + str(coef0))
    print("    Decision Function Shape = " + str(decision_function_shape))
    print("    Tolerance = " + str(tolerance))
    print("    Maximum Iteration = " + str(int(max_iteration)))
    print("    Breaking Prediction ties (break_ties) = " + str(break_ties))   
    print("    Classifier File = " + str(oclass))
    print("    Classifier File Type = " + str(otype))
    print("")
    classifier = SVC(kernel=str(kernel), 
                     degree=int(degree),
                     gamma=str(gamma),
                     coef0=float(coef0),
                     decision_function_shape=str(decision_function_shape),
                     tol=float(tolerance),
                     max_iter=int(max_iteration),
                     break_ties=break_ties)
    process_classifier(datafile, label, classifier, oclass, otype, 
                       classparam, confusion, classreport, cross_validation)
    print("===================== SVM Generated =====================")

def useSVM(datafile, classfile, classtype, resultfile):
    """!
    Function to use a previously generated support vector machine (SVM) 
    to classify data.

    Usage:

        python bactclass.py useSVM --datafile=classifier_use.csv --classfile=classifier_SVM.pickle --classtype=pickle --resultfile=classifier_result.csv
    
    @param datafile String: Path to CSV file containing data to be classified.
    @param classfile String: Path to the generated classifier.
    @param classtype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib".
    @param resultfile String: Path to write out the classified results.
    """
    useScikitClassifier("SVM", datafile, classfile, classtype, resultfile)

def generateDT(datafile, label, 
               oclass="classifier_DT.pickle", 
               otype="pickle",
               criterion="gini",
               splitter="best",
               max_depth=0,
               min_samples_split=2,
               min_samples_leaf=1,
               min_weight_fraction_leaf=0.0,
               max_leaf_nodes=0,
               ccp_alpha=0.0,
               classparam=True, 
               confusion=True, 
               classreport=True,
               cross_validation=5):
    """!
    Function to generate a decision tree from given data. 

    Usage:
        
        python bactclass.py genDT --datafile=classifier_train.csv --label=Class --oclass=classifier_DT.pickle --otype=pickle --criterion=gini --splitter=best --max_depth=0 --min_samples_split=2 --min_samples_leaf=1 --min_weight_fraction_leaf=0.0 --max_leaf_nodes=0 --ccp_alpha=0.0 --classparam=True --confusion=True --classreport=True --cross_validation=5

    @param datafile String: Path to CSV data file used to generate SVM.
    @param label String: Column (field) name in the data file to indicate the class label.
    @param oclass String: Path to write out the generated classifier. Default = classifier_SVM.pickle
    @param otype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib". Default = pickle
    @param criterion String: The function to measure the quality of a split. Allowable types are "gini" (Gini impurity) and "entropy" (information gain). Default = gini
    @param splitter String: The strategy used to choose the split at each node. Allowable types are "best" (best split) and "random" (best random split). Default = best
    @param max_depth Integer: The maximum depth of the tree. If less than 1, then nodes are expanded until all leaves are pure or until all leaves contain less than min_samples_split samples. Default = 0
    @param min_samples_split Integer: The minimum number of samples required to split an internal node. Default = 2
    @param min_samples_leaf Integer: The minimum number of samples required to be at a leaf node. A split point at any depth will only be considered if it leaves at least min_samples_leaf training samples in each of the left and right branches. This may have the effect of smoothing the model, especially in regression. Default = 1
    @param min_weight_fraction_leaf Float: The minimum weighted fraction of the sum total of weights (of all the input samples) required to be at a leaf node. Default = 0.0
    @param max_leaf_nodes Integer: Grow a tree with max_leaf_nodes in best-first fashion. Best nodes are defined as relative reduction in impurity. If less than 1, then unlimited number of leaf nodes. Default = 0
    @param ccp_alpha Float: Complexity parameter used for Minimal Cost-Complexity Pruning. The subtree with the largest cost complexity that is smaller than ccp_alpha will be chosen. If zero, no pruning is performed. Default = 0.0
    @param classparam Boolean: Flag to indicate whether to print out SVM parameters. Default = True
    @param confusion Boolean: Flag to indicate whether to print out confusion matrix. Default = True
    @param classreport Boolean: Flag to indicate whether to print out classification report. Default = True
    @param cross_validation Integer: Number of cross validation (if any) to perform. If less than 2, cross validation will not be carried out. Default = 5
    """
    from sklearn.tree import DecisionTreeClassifier
    print("")
    print("Task: Generate Decision Tree (DT) Classifier")
    print("Parameters:")
    print("    Data File = " + str(datafile))
    print("    Classification Label = " + str(label))
    print("    Citerion = " + str(criterion))
    print("    Splitting Strategy = " + str(splitter))
    print("    Maximum Tree Depth = " + str(max_depth))
    print("    Minimum Samples to Split = " + str(min_samples_split))
    print("    Minimum Samples per Leaf = " + str(min_samples_leaf))
    print("    Minimum Weight Fraction per Leaf = " + str(min_weight_fraction_leaf))
    print("    Maximum number of Leaf Nodes = " + str(max_leaf_nodes))
    print("    Minimal Cost-Complexity Pruning Parameter = " + str(ccp_alpha))
    print("    Classifier File = " + str(oclass))
    print("    Classifier File Type = " + str(otype))
    print("")
    max_depth = int(max_depth)
    if max_depth == 0: max_depth = None
    max_leaf_nodes = int(max_leaf_nodes)
    if max_leaf_nodes == 0: max_leaf_nodes = None
    classifier = DecisionTreeClassifier(criterion=str(criterion),
                                        splitter=str(splitter),
                                        max_depth=max_depth,
                                        min_samples_split=int(min_samples_split),
                                        min_samples_leaf=int(min_samples_leaf),
                                        min_weight_fraction_leaf=float(min_weight_fraction_leaf),
                                        max_leaf_nodes=max_leaf_nodes,
                                        ccp_alpha=float(ccp_alpha))
    process_classifier(datafile, label, classifier, oclass, otype, 
                       classparam, confusion, classreport, cross_validation)
    print("===================== DT Generated =====================")
    

def useDT(datafile, classfile, classtype, resultfile):
    """!
    Function to use a previously generated decision tree (DT) to classify data.

    Usage:

        python bactclass.py useDT --datafile=classifier_use.csv --classfile=classifier_DT.pickle --classtype=pickle --resultfile=classifier_result.csv
    
    @param datafile String: Path to CSV file containing data to be classified.
    @param classfile String: Path to the generated classifier.
    @param classtype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib".
    @param resultfile String: Path to write out the classified results.
    """
    useScikitClassifier("DT", datafile, classfile, classtype, resultfile)
    
    
def generateBNB(datafile, label, 
                oclass="classifier_BNB.pickle", 
                otype="pickle",
                alpha=1.0,
                binarize=0.0,
                fit_prior=True,
                classparam=True, 
                confusion=True, 
                classreport=True,
                cross_validation=5):
    """!
    Function to generate a Bernoulli Naive Bayes classifier (BNB) from given data. 

    Usage:
        
        python bactclass.py genBNB --datafile=classifier_train.csv --label=Class --oclass=classifier_BNB.pickle --otype=pickle --alpha=1.0 --binarize=0.0 --fit_prior=True --classparam=True --confusion=True --classreport=True --cross_validation=5

    @param datafile String: Path to CSV data file used to generate SVM.
    @param label String: Column (field) name in the data file to indicate the class label.
    @param oclass String: Path to write out the generated classifier. Default = classifier_SVM.pickle
    @param otype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib". Default = pickle
    @param alpha Float: Additive (Laplace/Lidstone) smoothing parameter (0 for no smoothing). Default = 1.0
    @param binarize Float: Threshold for binarizing (mapping to booleans) of sample features. Default=0.0
    @param fit_prior Boolean: Flag to indicate whether to learn class prior probabilities or not. If False, a uniform prior will be used. Default = True
    @param classparam Boolean: Flag to indicate whether to print out SVM parameters. Default = True
    @param confusion Boolean: Flag to indicate whether to print out confusion matrix. Default = True
    @param classreport Boolean: Flag to indicate whether to print out classification report. Default = True
    @param cross_validation Integer: Number of cross validation (if any) to perform. If less than 2, cross validation will not be carried out. Default = 5
    """
    from sklearn.naive_bayes import BernoulliNB
    print("")
    print("Task: Generate Bernoulli Naive Bayes Classifier")
    print("Parameters:")
    print("    Data File = " + str(datafile))
    print("    Classification Label = " + str(label))
    print("    Additive Smoothing Parameter = " + str(alpha))
    print("    Binarizing Threshold = " + str(binarize))
    print("    Learn Prior Probabilities = " + str(fit_prior))
    print("    Classifier File = " + str(oclass))
    print("    Classifier File Type = " + str(otype))
    print("")
    classifier = BernoulliNB(alpha=float(alpha),
                             binarize=float(binarize),
                             fit_prior=fit_prior)
    process_classifier(datafile, label, classifier, oclass, otype, 
                       classparam, confusion, classreport, cross_validation)
    print("===================== BernoulliNB Generated =====================")
    
def useBNB(datafile, classfile, classtype, resultfile):
    """!
    Function to use a previously generated Bernoulli Naive Bayes (BNB) classifier to classify data.

    Usage:
        python bactclass.py useBNB --datafile=classifier_use.csv --classfile=classifier_BNB.pickle --classtype=pickle --resultfile=classifier_result.csv
    
    @param datafile String: Path to CSV file containing data to be classified.
    @param classfile String: Path to the generated classifier.
    @param classtype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib".
    @param resultfile String: Path to write out the classified results.
    """
    useScikitClassifier("BNB", datafile, classfile, classtype, resultfile)
    
def generateCNB(datafile, label, 
                oclass="classifier_CNB.pickle", 
                otype="pickle",
                alpha=1.0,
                fit_prior=True,
                classparam=True, 
                confusion=True, 
                classreport=True,
                cross_validation=5):
    """!
    Function to generate a Complement Naive Bayes classifier (CNB) from given data. 

    Usage:
        
        python bactclass.py genCNB --datafile=classifier_train.csv --label=Class --oclass=classifier_CNB.pickle --otype=pickle --alpha=1.0 --fit_prior=True --classparam=True --confusion=True --classreport=True --cross_validation=5

    @param datafile String: Path to CSV data file used to generate SVM.
    @param label String: Column (field) name in the data file to indicate the class label.
    @param oclass String: Path to write out the generated classifier. Default = classifier_CNB.pickle
    @param otype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib". Default = pickle
    @param alpha Float: Additive (Laplace/Lidstone) smoothing parameter (0 for no smoothing). Default = 1.0
    @param fit_prior Boolean: Flag to indicate whether to learn class prior probabilities or not. If False, a uniform prior will be used. Default = True
    @param classparam Boolean: Flag to indicate whether to print out SVM parameters. Default = True
    @param confusion Boolean: Flag to indicate whether to print out confusion matrix. Default = True
    @param classreport Boolean: Flag to indicate whether to print out classification report. Default = True
    @param cross_validation Integer: Number of cross validation (if any) to perform. If less than 2, cross validation will not be carried out. Default = 5
    """
    from sklearn.naive_bayes import ComplementNB
    print("")
    print("Task: Generate Complement Naive Bayes Classifier")
    print("Parameters:")
    print("    Data File = " + str(datafile))
    print("    Classification Label = " + str(label))
    print("    Additive Smoothing Parameter = " + str(alpha))
    print("    Learn Prior Probabilities = " + str(fit_prior))
    print("    Classifier File = " + str(oclass))
    print("    Classifier File Type = " + str(otype))
    print("")
    classifier = ComplementNB(alpha=float(alpha),
                              fit_prior=fit_prior)
    process_classifier(datafile, label, classifier, oclass, otype, 
                       classparam, confusion, classreport, cross_validation)
    print("===================== MultinomialNB Generated =====================")
    
def useCNB(datafile, classfile, classtype, resultfile):
    """!
    Function to use a previously generated Complement Naive Bayes (CNB) classifier to classify data.
    
    Usage:
        python bactclass.py useCNB --datafile=classifier_use.csv --classfile=classifier_CNB.pickle --classtype=pickle --resultfile=classifier_result.csv
    
    @param datafile String: Path to CSV file containing data to be classified.
    @param classfile String: Path to the generated classifier.
    @param classtype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib".
    @param resultfile String: Path to write out the classified results.
    """
    useScikitClassifier("CNB", datafile, classfile, classtype, resultfile)

def generateMNB(datafile, label, 
                oclass="classifier_MNB.pickle", 
                otype="pickle",
                alpha=1.0,
                fit_prior=True,
                classparam=True, 
                confusion=True, 
                classreport=True,
                cross_validation=5):
    """!
    Function to generate a Multinomial Naive Bayes classifier (MNB) from given data. 

    Usage:
        
        python bactclass.py genMNB --datafile=classifier_train.csv --label=Class --oclass=classifier_BNB.pickle --otype=pickle --alpha=1.0 --fit_prior=True --classparam=True --confusion=True --classreport=True --cross_validation=5

    @param datafile String: Path to CSV data file used to generate SVM.
    @param label String: Column (field) name in the data file to indicate the class label.
    @param oclass String: Path to write out the generated classifier. Default = classifier_SVM.pickle
    @param otype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib". Default = pickle
    @param alpha Float: Additive (Laplace/Lidstone) smoothing parameter (0 for no smoothing). Default = 1.0
    @param fit_prior Boolean: Flag to indicate whether to learn class prior probabilities or not. If False, a uniform prior will be used. Default = True
    @param classparam Boolean: Flag to indicate whether to print out SVM parameters. Default = True
    @param confusion Boolean: Flag to indicate whether to print out confusion matrix. Default = True
    @param classreport Boolean: Flag to indicate whether to print out classification report. Default = True
    @param cross_validation Integer: Number of cross validation (if any) to perform. If less than 2, cross validation will not be carried out. Default = 5
    """
    from sklearn.naive_bayes import MultinomialNB
    print("")
    print("Task: Generate Multinomial Naive Bayes Classifier")
    print("Parameters:")
    print("    Data File = " + str(datafile))
    print("    Classification Label = " + str(label))
    print("    Additive Smoothing Parameter = " + str(alpha))
    print("    Learn Prior Probabilities = " + str(fit_prior))
    print("    Classifier File = " + str(oclass))
    print("    Classifier File Type = " + str(otype))
    print("")
    classifier = MultinomialNB(alpha=float(alpha),
                               fit_prior=fit_prior)
    process_classifier(datafile, label, classifier, oclass, otype, 
                       classparam, confusion, classreport, cross_validation)
    print("===================== MultinomialNB Generated =====================")
    
def useMNB(datafile, classfile, classtype, resultfile):
    """!
    Function to use a previously generated Multinomial Naive Bayes classifier to classify data.
    
    Usage:
        python bactclass.py useMNB --datafile=classifier_use.csv --classfile=classifier_MNB.pickle --classtype=pickle --resultfile=classifier_result.csv
    
    @param datafile String: Path to CSV file containing data to be classified.
    @param classfile String: Path to the generated classifier.
    @param classtype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib".
    @param resultfile String: Path to write out the classified results.
    """
    useScikitClassifier("MNB", datafile, classfile, classtype, resultfile)

def generateGNB(datafile, label, 
                oclass="classifier_GNB.pickle", 
                otype="pickle",
                var_smoothing=1e-9,
                classparam=True, 
                confusion=True, 
                classreport=True,
                cross_validation=5):
    """!
    Function to generate a Gaussian Naive Bayes classifier (GNB) from given data. 

    Usage:
        
        python bactclass.py genGNB --datafile=classifier_train.csv --label=Class --oclass=classifier_GNB.pickle --otype=pickle --var_smoothing=1e-9 --classparam=True --confusion=True --classreport=True --cross_validation=5

    @param datafile String: Path to CSV data file used to generate SVM.
    @param label String: Column (field) name in the data file to indicate the class label.
    @param oclass String: Path to write out the generated classifier. Default = classifier_SVM.pickle
    @param otype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib". Default = pickle
    @param var_smoothing Float: Portion of the largest variance of all features that is added to variances for calculation stability. Default = 1e-9
    @param classparam Boolean: Flag to indicate whether to print out SVM parameters. Default = True
    @param confusion Boolean: Flag to indicate whether to print out confusion matrix. Default = True
    @param classreport Boolean: Flag to indicate whether to print out classification report. Default = True
    @param cross_validation Integer: Number of cross validation (if any) to perform. If less than 2, cross validation will not be carried out. Default = 5
    """
    from sklearn.naive_bayes import GaussianNB
    print("")
    print("Task: Generate Gaussian Naive Bayes Classifier")
    print("Parameters:")
    print("    Data File = " + str(datafile))
    print("    Classification Label = " + str(label))   
    print("    Variance Smoothing = " + str(var_smoothing))
    print("    Classifier File = " + str(oclass))
    print("    Classifier File Type = " + str(otype))
    print("")
    classifier = GaussianNB(var_smoothing=float(var_smoothing))
    process_classifier(datafile, label, classifier, oclass, otype, 
                       classparam, confusion, classreport, cross_validation)
    print("===================== GaussianNB Generated =====================")
    
def useGNB(datafile, classfile, classtype, resultfile):
    """!
    Function to use a previously generated Gaussian Naive Bayes classifier (GNB) classifier to classify data.

    Usage:
        python bactclass.py useGNB --datafile=classifier_use.csv --classfile=classifier_GNB.pickle --classtype=pickle --resultfile=classifier_result.csv
    
    @param datafile String: Path to CSV file containing data to be classified.
    @param classfile String: Path to the generated classifier.
    @param classtype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib".
    @param resultfile String: Path to write out the classified results.
    """
    useScikitClassifier("GNB", datafile, classfile, classtype, resultfile)


if __name__ == "__main__":
    exposed_functions = {"genANN": generateANN,
                         "genBNB": generateBNB,
                         "genCNB": generateCNB,
                         "genDT": generateDT,
                         "genGNB": generateGNB,
                         "genMNB": generateMNB,
                         "genSVM": generateSVM,
                         "recycle": recycle,
                         "useANN": useANN,
                         "useBNB": useBNB,
                         "useCNB": useCNB,
                         "useDT": useDT,
                         "useGNB": useGNB,
                         "useMNB": useMNB,
                         "useSVM": useSVM}
    fire.Fire(exposed_functions)
