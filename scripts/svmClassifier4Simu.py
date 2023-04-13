import pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
import joblib
import os


def readCNVs(sample_name):
    """Read a CNV list from bed file into a DataFrame and add labels
    - sample_name: such as AK1, HG002 etc.
    return: two DataFrames, one is true-positive CNVs, the other is false-positive CNVs
    """
    toolScore = {'smoove': 5, 'delly': 4, 'cnvkit': 3, 'cnvpytor': 2, 'mops': 1}

    cnv_tp = pd.read_csv("/data3/jhsun/project/CNVPipe/analysis-CNVSimulator/evaluation-svm/"+sample_name+".observTP.bed", sep='\t', names=['accumScore', 'depthScore', 'tools', 'toolNum', 'cnvfilter', 'goodScore'], usecols=[5,6,7,8,9,10])
    cnv_tp['label'] = ['T' for _ in range(len(cnv_tp))]
    cnv_tp['cnvfilter'] = [0 if x==True else 1 for x in cnv_tp.cnvfilter]
    cnv_tp['tools'] = [sum(toolScore[t] for t in set(x.split(','))) for x in cnv_tp.tools]

    cnv_fp = pd.read_csv("/data3/jhsun/project/CNVPipe/analysis-CNVSimulator/evaluation-svm/"+sample_name+".FP.bed", sep='\t', names=['accumScore', 'depthScore', 'tools', 'toolNum', 'cnvfilter', 'goodScore'], usecols=[5,6,7,8,9,10])
    cnv_fp['label'] = ['F' for _ in range(len(cnv_fp))]
    cnv_fp['cnvfilter'] = [0 if x==True else 1 for x in cnv_fp.cnvfilter]
    cnv_fp['tools'] = [sum(toolScore[t] for t in set(x.split(','))) for x in cnv_fp.tools]

    cnv = pd.concat([cnv_tp, cnv_fp])
    return cnv

def buildSVMModel(model_file):

    cnvs = pd.concat([readCNVs('sample'+str(i)) for i in range(31,37)])
    # cnvs = pd.concat([readCNVs('sample43'), readCNVs('sample44'), readCNVs('sample45'), readCNVs('sample46'), readCNVs('sample47'), readCNVs('sample48')])

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

    training_vectors = ['accumScore', 'depthScore', 'tools', 'toolNum', 'cnvfilter', 'goodScore']

    # train model and store
    # 200k: 30x - 0.9363; 10x - 0.9372; 5x - 0.9250; 1x - 0.9375; 0.5x - 0.9545; 0.1x - 0.9565
    # 50k: 30x - 0.9293; 10x - 0.9348; 5x - 0.9082; 1x - 0.9294; 0.5x - 1; 0.1x - 0.8750
    # 10k: 30x - 0.9403; 10x - 0.9204; 5x - 0.9036; 1x - 0.7368; 0.5x - 0.5556; 0.1x - 1.0000
    model_file = "/data3/jhsun/github-repo/CNVPipe/resources/SVM/cnv_svm_classifier_simu_200k_5x.pkl"
    if not os.path.exists(model_file):
        buildSVMModel(model_file)

    # predict and evaluate
    clf = joblib.load(model_file)

    sample13 = readCNVs('sample31')
    test_x = sample13[training_vectors]
    test_y = sample13['label']
    predictions = clf.predict(test_x)
    
    accuScore = accuracy_score(test_y, predictions)
    print(f"Accuracy score is {accuScore:.4f}")