#The following code takes and merges multiple text files
#1.Merges the training and the test sets to create one data set.
#2.Extracts only the measurements on the mean and standard deviation for each measurement. 
#3.Uses descriptive activity names to name the activities in the data set
#4.Appropriately labels the data set with descriptive variable names. 
#5.independent tidy data set with the average of each variable for each activity and each subject.

install.packages("data.table")
library(data.table)
setwd("~/Documents/Learn/Getting_Cleaning Data/UCI HAR Dataset") #set WD
#Merge the test and train set to get one large data set 

#read in train data
subject_train<-read.table("train/subject_train.txt") #subject ID
train1<-read.table("train/X_train.txt")
train2<-read.table("train/Y_train.txt") #activity ID number, to be matched to act_labels

#read in test data
subject_test<-read.table("test/subject_test.txt") #subject ID
test1<-read.table("test/X_test.txt")
test2<-read.table("test/Y_test.txt") #activity ID number, to be matched to act_labels

#read in features and labels
features<-read.table("features.txt",stringsAsFactors=F)[[2]]
act_labels<-read.table("activity_labels.txt",stringsAsFactors=F)

tt<-rbind(train1,test1)#merge train and test data 
colnames(tt)<-features #label columns with features file
##################################################################

#Extracts only the measurements on the mean and standard deviation
tt <- tt[, grep("mean|std", names(tt))] #grep() to identify columns following mean and std pattern

##################################################################
#Use descriptive names to name activity in the dataset
all_subjects<-rbind(subject_train,subject_test) #complete build of entire data set 
all_act<-rbind(train2,test2)
colnames(all_subjects)<-c("Subject")
Activity<-act_labels$V2[all_act$V1] #Y_test.txt or Y_train.txt corresponds to activity_labels.txt
complete_tt<-cbind(all_subjects,Activity,tt)

##################################################################
#Use descriptibe names to name variables in the dataset 
VarNam<-variable.names(complete_tt) #Setting variable names as VarNam 
VarNam<-gsub("BodyBody","Body",VarNam) #removed redundancy but did not change the names dramatically as to not lose meaning. Descriptions too long to spell out "time" and "frequency"
VarNam<-gsub("[(][)]-X"," X",VarNam)
VarNam<-gsub("[(][)]-Y"," Y",VarNam)
VarNam<-gsub("[(][)]-Z"," Z",VarNam)
VarNam<-gsub("[(][)]"," ",VarNam)
names(complete_tt)<-VarNam
##################################################################

#tidy dataset with the average of each variable for each activity and each subject
colnames(complete_tt)[3:81]->id_variables #differentiate from Subject/Activity and other variables 
tidy1<-melt(complete_tt,id=c("Subject","Activity"), measure.vars=id_variables)
tidy2<-dcast(tidy1, Subject+Activity~ variable, mean)
##################################################################
write.table(tidy2, "tidy_data.txt", row.names=FALSE)

