library(neurohcp)
make_aws_call(path_to_file = "/", 
              bucket = "hcp-openaccess",region = "us-east-1", access_key = "", 
              secret_key = "",
              lifetime_minutes = 5, query = NULL, verb = "GET", sign = TRUE)
subject_data<- read.csv("subject_list_hcp_relational.csv",header = TRUE)
subjects <- subject_data$subject
scans <- subject_data$lr 

func_output_dir = 'func_rel/raw'
physio_output_dir = 'physio_rel/raw'
events_output_dir = 'events_rel'
dir.create(file.path(func_output_dir), recursive=TRUE, showWarnings = FALSE)
dir.create(file.path(physio_output_dir), recursive=TRUE, showWarnings = FALSE)
dir.create(file.path(events_output_dir), recursive=TRUE, showWarnings = FALSE)



# Loop through subjects and pull physio and func data
for (n in 1:length(subjects)) {
		i = subjects[n]
		s = scans[n]
		print(paste(n,': ',i))
		download_hcp_file(paste("HCP_1200/",i,"/MNINonLinear/Results/tfMRI_RELATIONAL_", s, "/tfMRI_RELATIONAL_", s, ".nii.gz",sep=""),
                    destfile = paste(func_output_dir,'/',i,"_", s, "_relational.nii.gz", sep = ""), error=FALSE)
        download_hcp_file(paste("HCP_1200/",i,"/MNINonLinear/Results/tfMRI_RELATIONAL_", s, "/tfMRI_RELATIONAL_", s, "_Physio_log.txt",sep=""),
                    destfile = paste(physio_output_dir,'/',i,"_", s, "_relational_physio.txt", sep = ""), error=FALSE)
        download_hcp_dir(paste("HCP_1200/",i,"/MNINonLinear/Results/tfMRI_RELATIONAL_", s, "/EVs",sep=""),
                    outdir = paste(events_output_dir,'/',i,"_", s, "_EV", sep = ""))
}

