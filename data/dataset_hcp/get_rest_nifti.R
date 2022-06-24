library(neurohcp)
make_aws_call(path_to_file = "/", 
              bucket = "hcp-openaccess",region = "us-east-1", access_key = "", 
              secret_key = "",
              lifetime_minutes = 5, query = NULL, verb = "GET", sign = TRUE)
subject_data<- read.csv("subject_list_hcp.csv",header = TRUE)
subjects <- subject_data$subject
scans <- subject_data$lr

func_output_dir = 'func/raw'
func_fix_output_dir = 'func_fix/raw'
physio_output_dir = 'physio/raw'
anat_output_dir = 'anat/raw_fs_seg'
dir.create(file.path(func_output_dir), recursive=TRUE, showWarnings = FALSE)
dir.create(file.path(func_fix_output_dir), recursive=TRUE, showWarnings = FALSE)
dir.create(file.path(physio_output_dir), recursive=TRUE, showWarnings = FALSE)
dir.create(file.path(anat_output_dir), recursive=TRUE, showWarnings = FALSE)


# Loop through subjects and pull physio, func and segmentation data
for (n in 1:length(subjects)) {
  i = subjects[n]
  s = scans[n]
  print(paste(n,': ',i))
  download_hcp_file(paste("HCP_1200/",i,"/MNINonLinear/Results/rfMRI_REST1_", s, "/rfMRI_REST1_", s, ".nii.gz",sep=""),
                    destfile = paste(func_output_dir,'/',i,"_", s, "1_rest.nii.gz", sep = ""), error=FALSE)
  download_hcp_file(paste("HCP_1200/",i,"/MNINonLinear/Results/rfMRI_REST1_", s, "/rfMRI_REST1_", s, "_hp2000_clean.nii.gz",sep=""), 
                    destfile = paste(func_fix_output_dir,'/',i,"_", s, "1_rest.nii.gz", sep = ""), error=FALSE)
  download_hcp_file(paste("HCP_1200/",i,"/MNINonLinear/Results/rfMRI_REST1_", s, "/rfMRI_REST1_", s, "_Physio_log.txt",sep=""),
                    destfile = paste(physio_output_dir,'/',i,"_", s, "1_physio.txt", sep = ""), error=FALSE)
  download_hcp_file(paste("HCP_1200/",i,"/MNINonLinear/aparc.a2009s+aseg.nii.gz",sep=""),
                    destfile = paste(anat_output_dir,'/',i,"_fsseg.nii.gz", sep = ""), error=FALSE)
}

