library(neurohcp)
make_aws_call(path_to_file = "/", 
              bucket = "hcp-openaccess",region = "us-east-1", access_key = "", 
              secret_key = "",
              lifetime_minutes = 5, query = NULL, verb = "GET", sign = TRUE)
subject_data<- read.csv("subject_list_hcp.csv",header = TRUE)
subjects <- subject_data$subject

func_output_dir = 'func/raw'
physio_output_dir = 'physio/raw'
dir.create(file.path(func_output_dir), recursive=TRUE, showWarnings = FALSE)
dir.create(file.path(physio_output_dir), recursive=TRUE, showWarnings = FALSE)

# Loop through subjects and pull physio and func data
for (n in 1:length(subjects)) {
  i = subjects[n]
  print(paste(n,': ',i))
  if (n %% 2 == 0) {
    download_hcp_file(paste("HCP_1200/",i,"/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL.nii.gz",sep=""), destfile = paste(func_output_dir,'/',i,"_RL1_rest.nii.gz", sep = ""), error=FALSE)
    download_hcp_file(paste("HCP_1200/",i,"/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_Physio_log.txt",sep=""), destfile = paste(physio_output_dir,'/',i,"_RL1_physio.txt", sep = ""), error=FALSE)
  }
  else {
    download_hcp_file(paste("HCP_1200/",i,"/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR.nii.gz",sep=""), destfile = paste(func_output_dir,'/',i,"_LR1_rest.nii.gz", sep = ""), error=FALSE)
    download_hcp_file(paste("HCP_1200/",i,"/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Physio_log.txt",sep=""), destfile = paste(physio_output_dir,'/',i,"_LR1_physio.txt", sep = ""), error=FALSE)
  }
}

