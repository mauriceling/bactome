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

def generateANN(datafile, label, 
                oclass="classifier_ANN.pickle", 
                otype="pickle",
                classparam=True, 
                confusion=True, 
                classreport=True):
    """!
    Function to generate an artificial neural network (multi-layer perceptron classifier) from given data. 

    Usage:
        
        python bactclass.py genANN --datafile=classifier_train.csv --label=Class --oclass=classifier_ANN.pickle --otype=pickle --classparam=True --confusion=True --classreport=True

    @param datafile String: Path to CSV data file used to generate SVM.
    @param label String: Column (field) name in the data file to indicate the class label.
    @param oclass String: Path to write out the generated classifier. Default = classifier_ANN.pickle
    @param otype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib". Default = pickle
    @param classparam Boolean: Flag to indicate whether to print out SVM parameters. Default = True
    @param confusion Boolean: Flag to indicate whether to print out confusion matrix. Default = True
    @param classreport Boolean: Flag to indicate whether to print out classification report. Default = True
    """
    from sklearn.neural_network import MLPClassifier
    print("")
    print("Task: Generate Artificial Neural Network (ANN) Multi-Layer Perceptron Classifier")
    print("Parameters:")
    print("    Data File = " + str(datafile))
    print("    Classification Label = " + str(label))
    print("    Classifier File = " + str(oclass))
    print("    Classifier File Type = " + str(otype))
    print("")
    X, Y = readData(datafile, True, label)
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size = 0.20)
    classifier = MLPClassifier(random_state=1, max_iter=300)
    classifier.fit(X_train, Y_train)
    Y_pred = classifier.predict(X_test)
    saveModel(oclass, otype, classifier)
    if classparam:
        showClassifierParameters(classifier)
    if confusion:
        showConfusionMatrix(Y_test, Y_pred)
    if classreport:
        showClassificationReport(Y_test, Y_pred)
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
    print("")
    print("Task: Classifying using Artificial Neural Network (ANN)")
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

def generateSVM(datafile, label, 
                oclass="classifier_SVM.pickle", 
                otype="pickle",
                kernel="linear", 
                degree=3,
                classparam=True, 
                confusion=True, 
                classreport=True):
    """!
    Function to generate a support vector machine from given data. 

    Usage:
        
        python bactclass.py genSVM --datafile=classifier_train.csv --label=Class --oclass=classifier_SVM.pickle --otype=pickle --kernel=linear --degree=3 --classparam=True --confusion=True --classreport=True

    @param datafile String: Path to CSV data file used to generate SVM.
    @param label String: Column (field) name in the data file to indicate the class label.
    @param oclass String: Path to write out the generated classifier. Default = classifier_SVM.pickle
    @param otype String: Type of file to write out the generated classifier. Allowable types are "pickle" and "joblib". Default = pickle
    @param kernel String: Type of SVM kernel. Acceptable values are linear, poly, rbf, sigmoid. Default = linear
    @param degree Integer: Polynomial degree, only used in polynomial (poly) kernel. Default = 3
    @param classparam Boolean: Flag to indicate whether to print out SVM parameters. Default = True
    @param confusion Boolean: Flag to indicate whether to print out confusion matrix. Default = True
    @param classreport Boolean: Flag to indicate whether to print out classification report. Default = True
    """
    from sklearn.svm import SVC
    print("")
    print("Task: Generate Support Vector Machine (SVM) Classifier")
    print("Parameters:")
    print("    Data File = " + str(datafile))
    print("    Classification Label = " + str(label))
    print("    SVM Kernel Type = " + str(kernel))
    print("    Polynomial Degree = " + str(int(degree)))
    print("    Classifier File = " + str(oclass))
    print("    Classifier File Type = " + str(otype))
    print("")
    X, Y = readData(datafile, True, label)
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size = 0.20)
    classifier = SVC(kernel=str(kernel), degree=int(degree))
    classifier.fit(X_train, Y_train)
    Y_pred = classifier.predict(X_test)
    saveModel(oclass, otype, classifier)
    if classparam:
        showClassifierParameters(classifier)
    if confusion:
        showConfusionMatrix(Y_test, Y_pred)
    if classreport:
        showClassificationReport(Y_test, Y_pred)
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
    print("")
    print("Task: Classifying using Support Vector Machine (SVM)")
    print("Parameters:")
    print("    Data File = " + str(datafile))
    print("    Classifier File = " + str(classfile))
    print("    Result File = " + str(resultfile))
    print("")
    classifier = loadModel(classfile, classtype)
    data = readData(datafile, False)
    data["Prediction"] = classifier.predict(data)
    data.to_csv(resultfile, index=False)


if __name__ == "__main__":
    exposed_functions = {"genANN": generateANN,
                         "genSVM": generateSVM,
                         "useANN": useANN,
                         "useSVM": useSVM}
    fire.Fire(exposed_functions)
