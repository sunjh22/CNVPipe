import pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import joblib
import os


def readCNVs(sample_name):
    """Read a CNV list from bed file into a DataFrame and add labels
    - sample_name: such as AK1, HG002 etc.
    return: two DataFrames, one is true-positive CNVs, the other is false-positive CNVs
    """
    toolScore = {'smoove': 5, 'delly': 4, 'cnvkit': 3, 'cnvpytor': 2, 'mops': 1}

    cnv_tp = pd.read_csv("/data3/jhsun/project/CNVPipe/analysis-CNVSimulator/evaluation-svm/"+sample_name+".observTP.bed", sep='\t', names=['accumScore', 'depthScore', 'tools', 'toolNum', 'cnvfilter'], usecols=[5,6,7,8,9])
    cnv_tp['label'] = ['T' for _ in range(len(cnv_tp))]
    cnv_tp['cnvfilter'] = [0 if x==True else 1 for x in cnv_tp.cnvfilter]
    cnv_tp['tools'] = [sum(toolScore[t] for t in set(x.split(','))) for x in cnv_tp.tools]

    cnv_fp = pd.read_csv("/data3/jhsun/project/CNVPipe/analysis-CNVSimulator/evaluation-svm/"+sample_name+".FP.bed", sep='\t', names=['accumScore', 'depthScore', 'tools', 'toolNum', 'cnvfilter'], usecols=[5,6,7,8,9])
    cnv_fp['label'] = ['F' for _ in range(len(cnv_fp))]
    cnv_fp['cnvfilter'] = [0 if x==True else 1 for x in cnv_fp.cnvfilter]
    cnv_fp['tools'] = [sum(toolScore[t] for t in set(x.split(','))) for x in cnv_fp.tools]

    cnv = pd.concat([cnv_tp, cnv_fp])
    return cnv

def testSVM(sample_name, cross_ratio):
    cnvs = readCNVs(sample_name)
    clf = make_pipeline(StandardScaler(), SVC())

    train, test = train_test_split(cnvs, test_size=cross_ratio)
    training_vectors = ['accumScore', 'depthScore', 'tools', 'toolNum', 'cnvfilter']
    train_x = train[training_vectors]
    test_x = test[training_vectors]

    train_y = train['label']
    test_y = test['label']

    clf.fit(train_x, train_y)

    # predict and evaluate
    predictions = clf.predict(test_x)
    accuScore = accuracy_score(test_y, predictions)
    return(accuScore)

def buildSVMModel(model_file):
    # sample7 = readCNVs('sample7')
    # sample8 = readCNVs('sample8')
    # sample9 = readCNVs('sample9')
    # sample10 = readCNVs('sample10')
    # sample11 = readCNVs('sample11')
    # sample12 = readCNVs('sample12')

    # cnvs = pd.concat([sample7, sample8, sample9, sample10, sample11, sample12])

    cnvs = pd.concat([readCNVs('sample13'), readCNVs('sample14'), readCNVs('sample15'), readCNVs('sample16'), readCNVs('sample17'), readCNVs('sample18')])

    # build SVC
    clf = make_pipeline(StandardScaler(), SVC())

    # train SVC
    train_x = cnvs[training_vectors]
    train_y = cnvs['label']
    clf.fit(train_x, train_y)

    # store the trained SVM model
    joblib.dump(clf, model_file)
    print("Successfully trained and stored SVM model")


if __name__ == '__main__':

    # samples = ['NA12878-1', 'NA12878-2', 'CHM13', 'AK1', 'HG002', 'HG00514', 'HG00733', 'NA19240']
    # sample = 'AK1'
    # cross_ratio = 0.3
    # for sample in samples:
    #     accuScore = testSVM(sample, cross_ratio=cross_ratio)
    #     print(f"Accuracy score for {sample} at cross_ratio {cross_ratio} is {accuScore:.4f}")

    # train, test = train_test_split(cnvs, test_size=0.2)
    # train_x = train[training_vectors]
    # test_x = test[training_vectors]
    # train_y = train['label']
    # test_y = test['label']

    # train_x = cnvs[training_vectors]
    # test_x = NA19240[training_vectors]
    # train_y = cnvs['label']
    # test_y = NA19240['label']
    # clf.fit(train_x, train_y)

    training_vectors = ['accumScore', 'depthScore', 'tools', 'toolNum', 'cnvfilter']

    # train model and store
    model_file = "/data3/jhsun/github-repo/CNVPipe/resources/SVM/cnv_svm_classifier_simu_30x.pkl"
    if not os.path.exists(model_file):
        buildSVMModel(model_file)

    # predict and evaluate
    clf = joblib.load(model_file)

    sample13 = readCNVs('sample14')
    test_x = sample13[training_vectors]
    test_y = sample13['label']
    predictions = clf.predict(test_x)
    
    accuScore = accuracy_score(test_y, predictions)
    print(f"Accuracy score is {accuScore:.4f}")