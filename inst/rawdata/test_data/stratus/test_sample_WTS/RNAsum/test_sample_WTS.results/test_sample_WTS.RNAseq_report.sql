### MySQL command to insert data for sample "test_sample_WTS"
use piedb;
INSERT INTO RNAseq_reports ( ID ,Platform, PatientID, SampleID, Cancer, Source, Project, Report, PMID, Analysis, Summary, Date ) VALUES ( 1000000, "RNA_seq", "", "test_sample_WTS", "TEST", "-", "-", "test_sample_WTS.RNAseq_report.html", "test_sample_WTS", "Findings summary,Fusion genes,Immune markers,HRD genes,Cancer genes,All genes,Drug matching,Input data", "Transcriptome summary for sample test_sample_WTS generated on 2022-01-21", "2022-01-21" )
  ON DUPLICATE KEY UPDATE ID=1000000 ,Platform="RNA_seq", PatientID="", SampleID="test_sample_WTS", Cancer="TEST", Source="-", Project="-", Report="test_sample_WTS.RNAseq_report.html", PMID="test_sample_WTS", Analysis="Findings summary,Fusion genes,Immune markers,HRD genes,Cancer genes,All genes,Drug matching,Input data", Summary="Transcriptome summary for sample test_sample_WTS generated on 2022-01-21", Date="2022-01-21";
SET @ID := 0;
UPDATE RNAseq_reports SET ID = ( SELECT @ID := @ID + 1 );
