
def data_mining():

    """This is a data mining project where we want to compared which machine learning technique is the best to classify the cell-communication (1) and non-cell-communication (0).

    The dataset is provided in DM_Project_24.csv with 106 columns and ~1600 rows.

    In here, we compared different machine learning techniques such as decision tree, random forest, naive bayes, k-nearest neighbour, adaboost, and voting classifier.
    """
    
    ### DATA PRE-PROCESSING ###

    #Import dataset
    import pandas as pd
    from sklearn import metrics
    from sklearn.model_selection import cross_val_score, KFold
    from sklearn.preprocessing import StandardScaler, MinMaxScaler
    from sklearn.impute import SimpleImputer
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import make_scorer, f1_score, accuracy_score, precision_score 

    """Don't forget to install panda and sklearn on the terminal."""

    # Reading the CSV
    data  = pd.read_csv('DM_Project_24.csv')

    # Statistical information 
    print(data.describe()) 
    #Lot of missing data found because total count only ~1400-1500

    # Splitting the datasets 
    X_train = data.iloc[:, :-1]  # Features: columns 0-105
    y_train = data.iloc[:, -1]   # Target: column 106

    # Columns: numerical and categorical
    numerical_cols = X_train.columns[:-2].tolist()  # First 103 columns
    categorical_cols = X_train.columns[-2:].tolist() # Last 2 columns

    # Open the test data 
    test_data = pd.read_csv("test_data.csv")

    # Using random forest as the model  
    rf = RandomForestClassifier(random_state=42)

    ### DATA PRE-PROCESSING ### 
    """More information can be seen in cleaning.py"""

    ## Imputation by class-specific ## 
    # Separate the data based on class
    class_groups = data.groupby(y_train)

    #Make list for new data 
    imputed_data = []
    for label, group in class_groups:
        #Imputation for numerical 
        numerical_imputer = SimpleImputer(strategy = "mean")
        group[numerical_cols] = numerical_imputer.fit_transform(group[numerical_cols])

        #Imputation for categorical 
        categorical_imputer = SimpleImputer(strategy='most_frequent')
        group[categorical_cols] = categorical_imputer.fit_transform(group[categorical_cols])
        
        imputed_data.append(group)
        
    #Make the list to a dataframe for future analysis 
    X_imputed_class = pd.concat(imputed_data).sort_index()
    print(X_imputed_class)

    #Cross-validation
    X_imputed_class = X_imputed_class.iloc[:, :-1]  # Features: columns 0-105
    y_imputed_class = X_imputed_class.iloc[:, -1]

    #Calculate F1
    cv_impu = cross_val_score(rf, X_imputed_class, y_imputed_class, cv=KFold(n_splits=5), scoring= "f1")
    print(f'F1 score across folds: {cv_impu.mean():.4f}')
    #F1 = 0.8217 for imputation 

    #Calculate accuracy 
    cv_impu_acc = cross_val_score(rf, X_imputed_class, y_imputed_class, cv=KFold(n_splits=5), scoring= "accuracy")
    print(f'Accuracy score across folds: {cv_impu_acc.mean():.4f}')
    #accuracy is 0.9621

    ## Normalization will used the dataset after imputation ##

    #Separate the numerical and categorical
    X = X_imputed_class.iloc[:, :-1]  # Features: columns 0-105
    y = X_imputed_class.iloc[:, -1]   # Target: column 106

    # Columns: numerical and categorical
    num_x_train = X.columns[:-2]  # First 103 columns
    nom_x_train = X.columns[-2:] # Last 2 columns

    # Separate numerical and categorical data
    X_num = X[num_x_train]
    X_cat = X[nom_x_train]

    # Ensure numerical columns are actually numeric
    X_num = X_num.apply(pd.to_numeric, errors='coerce')

    # Standardization (Z-score normalization)
    standard_scaler = StandardScaler()
    X_num_standardized = pd.DataFrame(
        standard_scaler.fit_transform(X_num),
        columns=X_num.columns,
        index=X_num.index)

    # Combine scaled numerical data with categorical data
    X_standardized = pd.concat([X_num_standardized, X_cat], axis=1)

    # See if it's succesfully normalized
    print(X_standardized)

    ### Cross validation combining the specific class and normalisation ###

    #Class specific imputation and z-score normalization
    #Calculate F1
    cv_standard_class = cross_val_score(rf, X_standardized, y_train, cv=KFold(n_splits=5), scoring= "f1")
    print(f'F1 score across folds: {cv_standard_class.mean():.4f}')
    #F1 = 0.7897  

    #Calculate accuracy 
    cv_acc_standard_c = cross_val_score(rf, X_standardized, y_train, cv=KFold(n_splits=5), scoring= "accuracy")
    print(f'Accuracy score across folds: {cv_acc_standard_c.mean():.4f}')
    #accuracy is 0.9556

    #Standardization for specific class 

    # Columns: numerical and categorical
    num_x_train = test_data.columns[:-2]  # First 103 columns
    nom_x_train = test_data.columns[-2:] # Last 2 columns

    # Separate numerical and ategorical data
    X_num = test_data[num_x_train]
    X_cat = test_data[nom_x_train]

    # Ensure numerical columns are actually numeric
    X_num = X_num.apply(pd.to_numeric, errors='coerce')

    # Standardization (Z-score normalization)
    standard_scaler = StandardScaler()
    X_num_standardized = pd.DataFrame(
        standard_scaler.fit_transform(X_num),
        columns=X_num.columns,
        index=X_num.index)

    # Make a new dataset
    test_standardized = pd.concat([X_num_standardized, X_cat], axis=1)
    test_standardized.to_csv('test_standard.csv', index=False)

    # Outliers detection using local outliers factor 

    # Importing extra library
    from sklearn.neighbors import LocalOutlierFactor

    lof = LocalOutlierFactor(n_neighbors=20, contamination=0.01)  # Adjust neighbors and contamination
    outliers = lof.fit_predict(X_standardized)
    mask = outliers == 1
    x_clean_lof = X_standardized[mask]
    y_clean_lof = y_train[mask]

    # Cross validation 
    # Accuracy
    cv_if_acc = cross_val_score(rf, x_clean_lof, y_clean_lof, cv=KFold(n_splits=5), scoring="accuracy")
    print(f'Accuracy score across folds: {cv_if_acc.mean():.3f}')
    #accuracy = 0.9621

    # F1
    cv_if = cross_val_score(rf, x_clean_lof, y_clean_lof, cv=KFold(n_splits=5), scoring= "f1")
    print(f'F1 score across folds: {cv_if.mean():.3f}')
    #F1 = 0.8217 

    #Final data have 1584 rows x 106 columns after removal of outliers

    #Combining the data for clean datasets -> if you to make a new datasets
    clean_data = pd.concat([x_clean_lof, y_clean_lof], axis=1)
    print(clean_data)
    clean_data.to_csv('train_data.csv', index=False)

    ### HYPERPARAMETER TUNING & MODEL SELECTION ###
    # Importing libraries for hyperparameter tuning  
    from sklearn.model_selection import GridSearchCV

    #Have a new clean dataset
    data = pd.read_csv("train_data.csv")

    # Splitting the datasets 
    X_train = data.iloc[:, :-1]  # Features: columns 0-105
    y_train = data.iloc[:, -1]   # Target: column 106

    # Columns: numerical and categorical
    num_cols = X_train.columns[:-2]  # First 103 columns
    nom_cols = X_train.columns[-2:] # Last 2 columns

    #Open the clean test data 
    test_data = pd.read_csv("test_standard.csv")

    ### Hyperparameter tuning for classifier ###

    ## Decision tree ##
    from sklearn.tree import DecisionTreeClassifier

    # Make the model
    dt = DecisionTreeClassifier(random_state =42)

    #Parameter grid 
    param_grid_dt = {
        'max_depth': [None, 10, 20, 30],
        'min_samples_split': [2, 5, 10],
        'min_samples_leaf': [1, 2, 4],
        'criterion': ['gini', 'entropy'],
        'max_features': ['sqrt', 'log2', None]}

    # Create the grid search object
    grid_search_dt = GridSearchCV(
        dt,
        param_grid_dt,
        cv=5,
        scoring="f1",
        n_jobs=-1,
        verbose=1)

    # Fit the grid search
    grid_search_dt.fit(X_train, y_train)

    # Print the best parameters and score
    print("Best parameters:", grid_search_dt.best_params_)
    print("Best cross-validation score:", grid_search_dt.best_score_)

    #Best parameters: {'criterion': 'gini', 'max_depth': 10, 'max_features': None,'min_samples_leaf': 4, 'min_samples_split': 2}
    #Best cross-validation score: 0.7958167645375827

    #Making predictions 
    df_tune = DecisionTreeClassifier(random_state=42, max_features=None, max_depth=10, 
                                    min_samples_leaf=4, min_samples_split=2)
    model_dt = df_tune.fit(X_train, y_train)

    #Evaluation or cross validation
    # F1 
    f1_dt = cross_val_score(model_dt, X_train, y_train, cv=5, scoring= "f1")
    print(f'The cross-validation f1-score is: {f1_dt.mean():.4f}')
    #The cross-validation f1-score is: 0.7958

    # Accuracy 
    acc_dt = cross_val_score(model_dt,  X_train, y_train, cv=5, scoring= 'accuracy')
    print(f'The accuracy: {acc_dt.mean():.4f}')
    #The accuracy: 0.9552

    ## Random forest ## 
    #Importing extra library
    from sklearn.ensemble import RandomForestClassifier

    #Classifier for random forest 
    rf = RandomForestClassifier(random_state = 42)

    #Grid CV 
    #Parameter grid 
    param_grid_rf ={ 
        'n_estimators': [100, 200, 300, 400, 500],
        'max_features': ['sqrt', 'log2',None],
        'max_depth' : [5, 7, 10],
        'criterion' :['gini', 'entropy']}

    # Create the grid search object
    grid_search_rf = GridSearchCV(
        rf,
        param_grid_rf,
        cv=5,
        scoring="f1",
        n_jobs=-1)

    # Fit the grid search
    grid_search_rf.fit(X_train, y_train)

    # Print the best parameters and score
    print("Best parameters:", grid_search_rf.best_params_)
    print("Best cross-validation score:", grid_search_rf.best_score_)

    #Best parameters: {'criterion': 'entropy', 'max_depth': 5, 'max_features': None, 'n_estimators': 100}
    #Best cross-validation score: 0.8355926251097454

    #Making predictions 
    rf_tune = RandomForestClassifier(random_state = 42, criterion='entropy', max_depth=5,
                                    max_features= None, n_estimators=100)

    rf_model = rf_tune.fit(X_train, y_train)

    #Evaluation training set
    #Calculation of f1
    f1_rf = cross_val_score(rf_model, X_train, y_train, cv=5, scoring= 'f1')
    print(f'The cross-validation f1-score is: {f1_rf.mean():.4f}')
    #The cross-validation f1-score is: 0.8356

    #Calculation of accuracy 
    acc_rf = cross_val_score(rf_model,  X_train, y_train, cv=5, scoring= 'accuracy')
    print(f'The accuracy: {acc_rf.mean():.4f}')
    #The accuracy: 0.9653

    ## Naive Bayes ## 
    from sklearn.naive_bayes import GaussianNB

    #Classifier or the model
    gnb = GaussianNB()

    #parameter grid 
    param_grid_nb = {'var_smoothing': [1e-9, 1e-8, 1e-7, 1e-6]}

    grid_search_nb= GridSearchCV(
        estimator=gnb, 
        param_grid=param_grid_nb, 
        scoring="f1", 
        cv=5)

    # Fit the grid search
    grid_search_nb.fit(X_train, y_train)

    # Print the best parameters and score
    print("Best parameters:", grid_search_nb.best_params_)
    print("Best cross-validation score:", grid_search_nb.best_score_)

    # Best parameters: {'var_smoothing': 1e-09}
    # Best cross-validation score: 0.2173875806245384

    ## K-Nearest Neighbour ##
    from sklearn.neighbors import KNeighborsClassifier

    #model 
    knn = KNeighborsClassifier()

    #Parameters
    param_grid_knn = {
        'n_neighbors': [3, 5, 7, 9],        # Number of neighbors to try  
        'metric': ['euclidean', 'manhattan', 'minkowski']  # Distance metrics
    }

    # Set up GridSearchCV with 5-fold cross-validation
    grid_search_knn = GridSearchCV(
        estimator=knn, 
        param_grid=param_grid_knn, 
        scoring='f1', 
        cv=5)

    # Fit the model
    grid_search_knn.fit(X_train, y_train)

    # Display the best parameters and accuracy score
    print("Best Parameters:", grid_search_knn.best_params_)
    print("Best Score:", grid_search_knn.best_score_)

    #Best Parameters: {'metric': 'manhattan', 'n_neighbors': 3}
    #Best Score: 0.2893607928714311

    ## Ensemble learning ## 
    ## Adaboost ##
    from sklearn.ensemble import AdaBoostClassifier
    ada = AdaBoostClassifier(random_state=42, algorithm='SAMME')

    #Parameters 
    param_grid_adaboost = {
        'n_estimators': [50, 100, 200],          
        'learning_rate': [0.1, 0.5, 1.0]}         

    #Grid search CV
    grid_search_ada = GridSearchCV(
            estimator=ada, 
            param_grid=param_grid_adaboost, 
            cv=5,
            scoring="f1")
        
    #fit the model 
    grid_search_ada.fit(X_train, y_train)

    # Display the best parameters and accuracy score
    print("Best Parameters:", grid_search_ada.best_params_)
    print("Best Score:", grid_search_ada.best_score_)

    #Best Parameters: {'learning_rate': 0.5, 'n_estimators': 200}
    #Best Score: 0.8373094035364523

    ada_tune = AdaBoostClassifier(random_state = 42, algorithm='SAMME', 
                                  learning_rate= 0.5, n_estimators= 200)

    model_ada = ada_tune.fit(X_train, y_train)

    #Evaluation 
    #F1 score
    f1_ada = cross_val_score(model_ada, X_train, y_train, cv=5, scoring= 'f1')
    print(f'The cross-validation f1-score is: {f1_ada.mean():.4f}')
    #The cross-validation f1-score is: 0.8373

    #Accuracy 
    acc_ada = cross_val_score(model_ada,  X_train, y_train, cv=5, scoring= 'accuracy')
    print(f'The accuracy: {acc_ada.mean():.4f}')
    #The accuracy: 0.9646

    ## Voting Classifier ##
    from sklearn.ensemble import VotingClassifier

    rf_model = RandomForestClassifier(n_estimators=100, random_state=42)
    dt_model = DecisionTreeClassifier(random_state=42)
    knn = KNeighborsClassifier()

    voting_clf = VotingClassifier(
        estimators=[('rf', rf_model), ('dt', dt_model), ('knn', knn)],
        voting='soft')

    # Define the parameter grid for GridSearchCV
    param_grid_vote = {
        'rf__n_estimators': [50, 100],
        'rf__max_depth': [5, 10],
        'dt__max_depth': [None, 10, 20, 30],
        'dt__min_samples_split': [2, 5, 10],
        'knn__n_neighbors': [3, 5, 7, 9],        # Number of neighbors to try  
        'knn__metric': ['euclidean', 'manhattan', 'minkowski'],
        'voting': ['soft', 'hard']
    }

    #Grid search 
    grid_search_vote = GridSearchCV(estimator=voting_clf, 
                           param_grid=param_grid_vote, 
                           cv=5, 
                           scoring="f1")
    #Fit to the model 
    grid_search_vote.fit(X_train, y_train)

    # Display the best parameters and accuracy score
    print("Best Parameters:", grid_search_vote.best_params_)
    print("Best Cross-Validation Score:", grid_search_vote.best_score_)

    #Best Parameters: {'dt__max_depth': 20, 'dt__min_samples_split': 2, 
    # 'knn__metric': 'manhattan', 'knn__n_neighbors': 9, 
    # 'rf__max_depth': 10, 'rf__n_estimators': 50, 'voting': 'soft'}
    #Best Cross-Validation Score: 0.8274940168610824

    #Fit the test data
    voting_clf.fit(X_train, y_train)

    #Cross validations
    #f1_scorer = make_scorer(f1_score, average='binary')
    f1_vote = cross_val_score(voting_clf, X_train, y_train, cv=5, scoring='f1')
    print(f'The cross-validation f1-score is: {f1_vote.mean():.4f}')

    #Accuracy 
    acc_vote = cross_val_score(voting_clf, X_train, y_train, cv=5, scoring='accuracy')
    print(f'The cross-validation acc-score is: {acc_vote.mean():.4f}')

    #The cross-validation f1-score is: 0.8093
    #The cross-validation acc-score is: 0.9577

    ### PREDICTION ###
    # The chosen model is adaboost due to the highest CV 
    model_ada = ada_tune.fit(X_train, y_train)
    y_pred_ada = ada_tune.predict(test_data)

    #There are 86 samples scored one - classified has cell-communication 

    ### GENERATION OF RESULT REPORT ###
    result = []
    for truth in y_pred_ada:
        result.append(truth)
    
    #F1 score - round it
    mean_f1 = round(f1_ada.mean(), 3)
    #Accuracy score - round it 
    mean_acc = round(acc_ada.mean(), 3)
    
    #Add the result to the list 
    result.append(mean_acc)
    result.append(mean_f1)

    print(result)

def main():
    data_mining()

if __name__ == "__main__":
    main()
    


