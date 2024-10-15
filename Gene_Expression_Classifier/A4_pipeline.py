

#####

#import libraries/packages/modules
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_predict
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
import seaborn as sns

##########------------------- Part 1 --------------------###########
print('--- Part 1 ---')
##Step 1: read metadata file and platform file to obtain gene symbols
print('Reading the metadata and platform file...')

metadata = pd.read_csv('meta_data.csv')
platform_df = pd.read_csv('GPL570-55999.txt', sep='\t', index_col=0, header=16, low_memory=False)
platform_df.reset_index(inplace=True)

#filter metadata to extract sample IDs for flu and rsv samples
#using '.tolist()' function create respective lists containing the samples IDs for flu and rsv samples
print('Filtering sample types in metadata...')
flu_samples = metadata[metadata['infection_status'] == 'influenza']['Sample_geo_accession'].tolist()
rsv_samples = metadata[metadata['infection_status'] == 'rsv']['Sample_geo_accession'].tolist()
control_samples = metadata[metadata['infection_status'] == 'none']['Sample_geo_accession'].tolist()

##Step 2: read matrix csv file
print('Reading matrix file and creating dataframe...')

#df = pd.read_csv('GSE34205_series_matrix_clean.txt', sep='\t', nrows=1000) #looking at 1000 rows
df = pd.read_csv('GSE34205_series_matrix_clean.txt', sep='\t', index_col=0) #looking at entire dataset

#reorganise columns by sample type (flu, rsv, Control), by extracting the columns directly from the dataframe
print('Reorganising dataframe and adding columns...')
#reorder columns to group flu, rsv and Control samples together
all_samples = flu_samples+rsv_samples+control_samples #create a new order
df = df.reindex(columns=all_samples) #reorder matrix

##Step 3: adding columns to dataframe

df_len = len(df)
df = df.assign(mean_rsv_ratio=np.zeros(df_len))
df = df.assign(mean_flu_ratio=np.zeros(df_len))
df = df.assign(p_value_rsv=np.zeros(df_len))
df = df.assign(corrected_p_value_rsv=np.zeros(df_len))
df = df.assign(p_value_flu=np.zeros(df_len))
df = df.assign(corrected_p_value_flu=np.zeros(df_len))
df["Selected_flu_feature"] = False
df["Selected_rsv_feature"] = False
df["Selected_feature"] = False

##Step 4: calculating mean and P-values
print('Calculating mean and perform t-test...')

#identify the values for each sample, log and mean
flu_values = np.mean(np.log(df[flu_samples]), axis=1)
rsv_values = np.mean(np.log(df[rsv_samples]), axis=1)

#calculate the mean of the log for each sample and add to respective column in dataframe
df['mean_flu_ratio'] = flu_values
df['mean_rsv_ratio'] = rsv_values

#perform 1-sample ttest of the log of the sample value in each row
df['p_value_flu'] = stats.ttest_1samp(np.log(df[flu_samples]), 0.0, axis=1).pvalue #0.0 is null hypothesis
df['p_value_rsv'] = stats.ttest_1samp(np.log(df[rsv_samples]), 0.0, axis=1).pvalue

##Step 5: Correct p-values with Bonferroni Correction for multiple testing

#multiply the p-value by number of rows
df['corrected_p_value_flu'] = df['p_value_flu'] * df_len
df['corrected_p_value_rsv'] = df['p_value_rsv'] * df_len

##Step 6: Label differentially expressed genes

#call genes with corrected p-value less than threshold & mean absolute value greater than 1 and set to boolean 'True'
print('Identifying and labelling differentially expressed genes...')
df['Selected_flu_feature'] = (df['corrected_p_value_flu'] < 0.05) & (abs(df['mean_flu_ratio']) >1)
df['Selected_rsv_feature'] = (df['corrected_p_value_rsv'] < 0.05) & (abs(df['mean_rsv_ratio']) >1)

#if either flu or rsv samples pass the conditions, set selected feature to boolean 'True'
df['Selected_feature'] = df['Selected_flu_feature'] | df['Selected_rsv_feature']

##Step 7: Save dataframes to .csv files
print('Saving dataframes as .csv files...')

#save entired dataframe with additional columns to .csv file
df.reset_index(inplace=True)
df.to_csv('matrix_plus_stats.csv', index=False)

#extract only rows with selected_features = True save to a second .csv file
selected_features = df[df['Selected_feature']]
selected_features.to_csv('features.csv', index=False)

##Step 8: Generate volcano plots of differentially expressed genes for both rsv and flu samples
print('Plotting two volcano plots...')

#filter dataframe for selected features
selected_features_flu = df[df['Selected_flu_feature']]
selected_features_rsv = df[df['Selected_rsv_feature']]

#merge the selected features with the platform df to obtain the gene symbols
selected_features_flu = pd.merge(selected_features_flu, platform_df, left_on='ID_REF', right_on='ID', how='left')
selected_features_rsv = pd.merge(selected_features_rsv, platform_df, left_on='ID_REF', right_on='ID', how='left')

####
#Volcano of Flu samples#
#create subplots and axis
fig, ax = plt.subplots(figsize=(8, 12), dpi=200)
ax.set_facecolor('whitesmoke') #plot face colour
ax.spines[['top', 'bottom', 'left', 'right']].set_color('none')  #remove border

#plotting the mean and -log p-value of all flu-samples
ax.scatter(df['mean_flu_ratio'], -np.log(df['corrected_p_value_flu']),
           s=5, color='lightgreen', label='Non-significant', edgecolor='black', linewidths=0.1, alpha=0.5)

#plotting mean and -log p-value of differentially expressed genes from flu samples with gene symbols as labels
for i, row in selected_features_flu[selected_features_flu['Selected_flu_feature']].iterrows():
    ax.scatter(row['mean_flu_ratio'], -np.log(row['corrected_p_value_flu']),
               s=10, color='green', edgecolor='black', linewidths=0.6, alpha=0.7, label=None)  #no label for individual points
    ax.text(row['mean_flu_ratio'] +0.05, -np.log(row['corrected_p_value_flu']), row['Gene Symbol'], fontsize=3.5)

#let Matplotlib automatically determine y-axis limits
ax.autoscale(enable=True, axis='y')

ax.set_xlabel('Mean Log Ratio')  #x-axis label
ax.set_ylabel('-Log Corrected p-value')  #y-axis label
ax.set_title('Volcano Plot for Differentially Expressed Genes - Flu Samples')  #title
ax.legend(labels=['Non-significant', 'Significant'])  #display legend with specific labels
ax.grid(True, color='lightgrey', linestyle='--')  #add grid to figure and color grid lines

plt.tight_layout() #prevent overlapping
#save plot as .png, set resolution and adjust bounding box to include all points without extra whitespace
plt.savefig('volcano_flu.png', dpi=300, bbox_inches='tight')
plt.show()  # display plot on screen

####
#Volcano of RSV samples#
#create subplots and axis
fig, ax = plt.subplots(figsize=(8, 12), dpi=200)
ax.set_facecolor('whitesmoke')
ax.spines[['top', 'bottom', 'left', 'right']].set_color('none')  #remove border

#plotting the mean and -log p-value of all rsv-samples
ax.scatter(df['mean_rsv_ratio'], -np.log(df['corrected_p_value_rsv']),
           s=5, color='lavender', label='Non-significant', edgecolor='black', linewidths=0.1, alpha=0.5)

#plotting mean and -log p-value of differentially expressed genes from rsv samples with gene symbols as labels
for i, row in selected_features_rsv[selected_features_rsv['Selected_rsv_feature']].iterrows():
    ax.scatter(row['mean_rsv_ratio'], -np.log(row['corrected_p_value_rsv']),
               s=10, color='blue', edgecolor='black', linewidths=0.6, alpha=0.7, label=None)  # No label for individual points
    ax.text(row['mean_rsv_ratio'] +0.05, -np.log(row['corrected_p_value_rsv']), row['Gene Symbol'], fontsize=3.5)

#let Matplotlib automatically determine y-axis limits
ax.autoscale(enable=True, axis='y')

ax.set_xlabel('Mean Log Ratio')  #x-axis label
ax.set_ylabel('-Log Corrected p-value')  #y-axis label
ax.set_title('Volcano Plot for Differentially Expressed Genes - RSV Samples')  #title
ax.legend(labels=['Non-significant', 'Significant'])  #display legend with specific labels
ax.grid(True, color='lightgrey', linestyle='--')  #add grid to figure and color grid lines

plt.tight_layout() #prevent overlapping
#save plot as .png, set resolution and adjust bounding box to include all points without extra whitespace
plt.savefig('volcano_rsv.png', dpi=300, bbox_inches='tight')
plt.show()  # display plot on screen



##########------------------- Part 2 --------------------###########
print('--- Part 2 ---')
##Step 1: Read in features.csv
print('Reading features.csv...')

features_df = pd.read_csv('features.csv', index_col=0)

##Step 2: Transpose the dataframe, making samples as rows and genes as columns with '.T'
print('Transposing and indexing features dataframe...')

transposed_df = features_df.iloc[:, 0:101].T #leave out the p-value, mean and Boolean columns
#add index 'Sample_geo_accession' into the dataframe
transposed_df.index = transposed_df.index.rename("Sample_geo_accession")

##Step 3: Standardise data to ensure zero mean and unit variance & generate 2 component PCA
#standardise data
print('Standardising data, performing PCA')

standardised_data = StandardScaler().fit_transform(transposed_df)

#apply PCA with 2 components
pca = PCA(n_components=2)
principal_components = pca.fit_transform(standardised_data)

#explained variance for each component
print(f'Explained Variance of PC1: {pca.explained_variance_ratio_[0]:.4f}')
print(f'Explained Variance of PC2: {pca.explained_variance_ratio_[1]:.4f}')

#total explained variance
total_variance = sum(pca.explained_variance_ratio_)
print(f'Total Explained Variance (PC1 + PC2): {total_variance:.4f}')

##Step 4: Create a new dataframe containing the PCA information
print('Generating new dataframe for PCA and merging with metadata...')

pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
pca_df.index = transposed_df.index #set index to match transposed df index for plotting

##Step 5: Merge the PCA dataframe with the metadata

pca_data = pd.merge(pca_df, metadata, left_index=True, right_on='Sample_geo_accession')

##Step 6: Plot the PCA plots
print('Plotting three PCA plots...')

#infection status PCA#
#create infection status PCA figure and axis
fig, ax = plt.subplots(figsize=(10, 8), dpi=200)
#set x & y-axis labels
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.grid(True, color='lightgrey', linestyle='--') #add grid to figure and colour grid lines

#define infection categories
infection_status = ['influenza', 'rsv', 'none'] #list for category labels
colours_infection = ['firebrick', 'seagreen', 'cornflowerblue'] #corresponding colours to category labels

#plot points based on infection status
for i in range(0,3): #3 infection status'
    temp_infection = infection_status[i] #current infection status
    colour_infection = colours_infection[i] #corresponding colour
    #new dataframe for specified infection status
    new_df_infection = pca_data[pca_data['infection_status'] == temp_infection]

    # scatterplot for current infection status
    ax.scatter(new_df_infection.loc[:, 'PC1'],
               new_df_infection.loc[:, 'PC2'],
               c=colour_infection,
               s=80,
               facecolors='darkgrey')

ax.legend(infection_status) #legend with infection status labels
fig.savefig('infection_pca.png', dpi=300) #save plot
plt.show() #display plot

####
#gender PCA#
#create gender PCA figure and axis
fig, ax = plt.subplots(figsize=(10, 8), dpi=200)
#set x & y-axis labels
ax.set_xlabel('Principal Component 1', fontsize=15)
ax.set_ylabel('Principal Component 2', fontsize=15)
ax.grid(True, color='lightgrey', linestyle='--') #add grid to plot

#define gender categories
genders = ['M', 'F'] #list for category labels
colours_gender = ['purple', 'orange'] #corresponding colours to category labels

#plot points based on gender
for i in range(0, 2): #2 genders
    temp_gender = genders[i] #current gender
    colour_gender = colours_gender[i] #corresponding colour
    #new dataframe for specified gender
    new_df_gender = pca_data[pca_data['gender'] == temp_gender]

    #scatterplot for current gender
    ax.scatter(new_df_gender.loc[:, 'PC1'],
                      new_df_gender.loc[:, 'PC2'],
                      c=colour_gender,
                      s=80,
                      facecolors='gray')

ax.legend(genders) #legend with gender group labels
fig.savefig('gender_pca.png', dpi=300) #save plot
plt.show()

####
#age group PCA#
#create age group PCA figure and axis
fig, ax = plt.subplots(figsize=(10, 8), dpi=200)
#set x & y-axis labels
ax.set_xlabel('Principal Component 1', fontsize=15)
ax.set_ylabel('Principal Component 2', fontsize=15)
ax.grid(True, color='lightgrey', linestyle='--') #add grid to plot

#define age groups
threshold = 6 #threshold for categorising
age_groups = ['> 6 months', '< 6 months'] #list for category labels
colours_age = ['brown', 'teal'] #corresponding colours to category labels

#plot points based on age group
for i in range(0, 2): #2 age groups
    temp_age = age_groups[i] #current age
    colour_age = colours_age[i] #corresponding colour
    #new dataframe for specified age group
    new_df_age = pca_data[pca_data['age_months'] == temp_age]

    #filter data based on specified age group
    if temp_age == '> 6 months':
        new_df_age = pca_data[pca_data['age_months'] >= threshold]
    elif temp_age == '< 6 months':
        new_df_age = pca_data[pca_data['age_months'] < threshold]

    #scatterplot for current age group
    ax.scatter(new_df_age.loc[:, 'PC1'],
                   new_df_age.loc[:, 'PC2'],
                   c=colour_age,
                   s=80,
                   facecolors='gray')

ax.legend(age_groups) #legend with age group labels
fig.savefig('age_group_pca.png', dpi=300) #save plot
plt.show()



##########------------------- Part 3 --------------------###########
print('--- Part 3 ---')
##Step 1: Prepare data for the machine learning
print('Preparing data for the machine learning...')

#features (x) is the transposed data
x = transposed_df
#labels (y) is the infection status from metadata
y = ['influenza'] * len(flu_samples) + ['rsv'] * len(rsv_samples) + ['none'] * len(control_samples)

##Step 2: Split data into training & test (70:30)

x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.3)

##Step 3: Create the Random Forest Classifer
print('Creating and training the Random Forest Classifier...')

classifier = RandomForestClassifier(n_estimators=100)

##Step 4: Train the model using the training data

classifier.fit(x_train,y_train)

##Step 5: Try model and make predictions on the test data
print('Testing classifier model...')

y_pred = classifier.predict(x_test)

##Step 6: Evaluate the model, create dataframes and save files
print('Evalutating the model, creating and saving reports...')

class_labels = ['influenza', 'RSV', 'none'] #list object of the sample types

#the confusion matrix - create dataframe and save as a .csv file
confuse_matrix = confusion_matrix(y_test, y_pred)
confuse_matrix_df = pd.DataFrame(confuse_matrix, index=class_labels, columns=class_labels)
confuse_matrix_df.to_csv('confusion_matrix.csv')

#create the classification report and save as a .txt file
class_report = classification_report(y_test, y_pred, target_names=class_labels)
with open('classification_report.txt', 'w') as report_file:
    report_file.write(class_report)

#the accuracy score
accuracy = accuracy_score(y_test, y_pred)

#feature importance - create dataframe and save as a .csv file
importance_df = pd.DataFrame({'ID': x.columns, 'Importance': classifier.feature_importances_})
importance_df = importance_df.sort_values(by='Importance', ascending=False) #sort to view the importance

#extract genbank IDs and gene symbols (for those genes that have an assigned symbol e.g. “MAPK1”)
#add these more meaningful gene symbols to the variable importance files
importance_df['GB_ACC'] = importance_df['ID'].map(platform_df.set_index('ID')['GB_ACC'])
importance_df['Gene Symbol'] = importance_df['ID'].map(platform_df.set_index('ID')['Gene Symbol'])

importance_df.to_csv('variable_importance.csv', index=False)

#print model evaluations
print("Confusion Matrix:\n", confuse_matrix)
print("\nClassification Report:\n", class_report)
print("\nAccuracy Score:", accuracy)

##Step 8: Visualising the initial classification
print('Creating a confusion matrix image...')

#confusion matrix
labels = ['influenza', 'rsv', 'none']
cm = confusion_matrix(y_test,y_pred,labels=labels)
fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(cm)
plt.title('Confusion matrix of the classifier')
fig.colorbar(cax)
ax.set_xticklabels([''] + labels) #these are misaligned unless you add an extra empty label
ax.set_yticklabels([''] + labels)
plt.xlabel('Predicted')
plt.ylabel('True')
plt.savefig("confused.png", dpi=300)

##Step 9: Run classifier multiple times and produce a summary report
print('Running classifier multiple times and producing a summary report...')

#create user input iterations object
#while true and try-except loop to ensure user inputs correct format for number of iterations
while True:
    try:
        iterations = int(input("Enter the number of iterations: "))
        if iterations <= 0:
            print("Please enter a positive integer.")
        else:
            break  #exit the loop if a valid input is provided
    except ValueError:
        print("Invalid input. Please enter a positive integer.")

results = [] #open list for the results
classification_reports = [] #create an empty list to store classification reports from each iteration

#perform cross-validated training and testing for the specified number of iterations
for _ in range(iterations):
    #train the classifier
    classifier.fit(x_train, y_train)

    #predict on the test set
    y_pred = classifier.predict(x_test)

    #collect information from each iteration
    iteration_results = {'Confusion_Matrix': confusion_matrix(y_test, y_pred),
                         'Classification_Report': classification_report(y_test, y_pred, target_names=class_labels, output_dict=True),
                         'Accuracy_Score': accuracy_score(y_test, y_pred)}

    #append information from each iteration into the results list previously opened
    results.append(iteration_results)

    #collect classification report for this iteration
    report_dict = classification_report(y_test, y_pred, target_names=class_labels, output_dict=True)

    #ensure that the classification report is a dictionary
    if isinstance(report_dict, dict):
        classification_reports.append(report_dict)
    else:
        print(f"Invalid classification report format.")


#calculate averages
avg_conf_matrix = np.mean([result['Confusion_Matrix'] for result in results], axis=0)

avg_accuracy_scores = np.mean([result['Accuracy_Score'] for result in results], axis=0)

#calculate the mean classification report
avg_class_report = {label: {
        metric: np.mean([result[label][metric] for result in classification_reports if isinstance(result[label][metric], (int, float))])
        for metric in ['precision', 'recall', 'f1-score', 'support']} for label in class_labels}

#write results to a summary file
with open('summary_results.txt', 'w') as summary_file:
    summary_file.write(f"Confusion Matrix:\n{avg_conf_matrix}\n")

#classification report formated with f-strings to create a tabular format in the file
#this is because the results are in a dictionary format currently
    #'header row'
    summary_file.write(f"{'':<15}{'precision':<15}{'recall':<15}{'f1-score':<15}{'support':<15}\n")

    #loop through each class label to write precision, recall, f1-score, and support to the file
    for label in class_labels: #class_labels are the sample types
        summary_file.write(
            f"{label:<15}" #class label to be a 'width of 15 characters'
            f"{avg_class_report[label]['precision']:.2f}{'':<15}"  #average precision formatted as a float with two decimal places
            f"{avg_class_report[label]['recall']:.2f}{'':<15}"  #average recall formatted as a float with two decimal places
            f"{avg_class_report[label]['f1-score']:.2f}{'':<15}"  #average f1-score formatted as a float with two decimal places
            f"{avg_class_report[label]['support']:.2f}\n")  #average support formatted as a float with two decimal places

    summary_file.write(f'Average Accuracy Score:\n{avg_accuracy_scores}\n')



##########------------------- Part 4 --------------------###########

#producing a large grid plot showing boxplots of the expression of the 100 most important features
print('Creating 100 boxplots...')
##Step 1: extract the top 100 most important features
top_feature_importance = importance_df.head(100)

##Step 2: merge transposed_df with metadata & create new boxplot dataframe

transposed_df.index=transposed_df.index.rename('Sample_geo_accession') #prep for merging
#merge to obtain the infection status
merged_df = pd.merge(transposed_df, metadata, left_index=True, right_on='Sample_geo_accession')
#create boxplot df with pd.melt to reshape in preparation for boxplot visualisation
boxplot_df = pd.melt(merged_df, id_vars=['infection_status'], value_vars=top_feature_importance['ID'],
                      var_name='Feature', value_name='Importance')

##Step 3: plot the boxplots and save as .png file

#plotting the boxplots, specifying 10*10 grid with shared axis
fig,axes=plt.subplots(10, 10, figsize=(20,20), sharex=True, sharey=True, dpi=200)
#flatten ax array for easier iteration with a single 'for in' loop
axes = axes.flatten()

#loop each feature (ID) through from the top features
for i, feature in enumerate(top_feature_importance['ID']):
    #create boxplot of current feature
    sns.boxplot(x='infection_status', y='Importance',
                data=boxplot_df[boxplot_df['Feature']==feature],
                ax=axes[i], showfliers=False) #exclude outliers
    axes[i].set_title(feature)
    axes[i].set_xlabel('') #remove individual x-axis labels
    axes[i].set_ylabel('') #remove individual y-axis labels

#set y-axis scale to logarithmic for all subplots
for ax in axes:
    ax.set_yscale('log')

plt.tight_layout() #prevent overlapping
plt.savefig('100_boxplots.png', dpi=300) #save boxplots
plt.show()

####
#producing a gene expression heatmap
print('Creating heatmaps...')

##Step 1: sort features_df to only consist of the gene expressions

#use only the relevant data i.e. gene expressions by dropping all the other columns
heatmap_data = features_df.drop(['mean_rsv_ratio', 'mean_flu_ratio', 'p_value_rsv', 'corrected_p_value_rsv', 'p_value_flu', 'corrected_p_value_flu',
                                   'Selected_flu_feature', 'Selected_rsv_feature', 'Selected_feature'], axis=1)

##Step 2: calculate the correlation matrix of the gene expression between samples

#using .corr()
heatmap_corr = heatmap_data.corr()

##Step 3: plot a heatmap and a clustermap then save both as .png files

#heatmap#
fig, ax = plt.subplots(figsize=(20, 15), dpi=200) #set figure size
#heatmap using the sns.heatmap function
sns_heatmap = sns.heatmap(heatmap_corr, cmap='coolwarm')
plt.savefig('heatmap.png', dpi=300) #save plot
plt.show()

#clustermap#
fig, ax = plt.subplots(figsize=(20, 15), dpi=200) #set figure size
#clustermap using the sns.clustermap function
sns_clustermap = sns.clustermap(heatmap_corr, cmap='coolwarm')
plt.savefig('clustermap.png', dpi=300) #save plot
plt.show()

print('Pipeline complete.')
