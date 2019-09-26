## Notes

1. Launch an instance using following information:
    - uisng `ami-0dc96254d5535925f`
    - `wts_report_ecsInstanceRole`
    - 100GB volume

2. Login to the instance: 

    - `wget https://github.com/umccr/RNAseq-Analysis-Report/blob/master/deploy/files/fetch-instance-bootstrap.sh` (place it in the _appropriate_ location).

3. Followed `Installing the Amazon ECS Container Agent on an Amazon Linux 2 EC2 Instance` from https://docs.aws.amazon.com/AmazonECS/latest/developerguide/ecs-agent-install.html
