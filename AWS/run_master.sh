cname=$1
compute_instance_type=$2
nNodes=$3
az=$4
tileList=$5

config=~/configs/${cname}

cp ~/.parallelcluster/config ${config}

#Set instance type and number of compute nodes
sed -i -e "s/.*compute_instance_type =.*/compute_instance_type = ${compute_instance_type}/" ${config}
sed -i -e "s/.*initial_queue_size =.*/initial_queue_size = ${nNodes}/" ${config}
sed -i -e "s/.*max_queue_size =.*/max_queue_size = ${nNodes}/" ${config}

#Set subnet in config file based on availability zone provided
if [ "${az}" = "a" ]; then subnet="subnet-827ff6ac";fi
if [ "${az}" = "b" ]; then subnet="subnet-92c597d8";fi
if [ "${az}" = "c" ]; then subnet="subnet-5f078803";fi
if [ "${az}" = "d" ]; then subnet="subnet-a95cdece";fi
if [ "${az}" = "f" ]; then subnet="subnet-bfb0f4b0";fi
sed -i -e "s/.*master_subnet_id =.*/master_subnet_id = ${subnet}/" ${config}




pcluster create -c ${config} ${cname}
ip=$(pcluster status $cname | grep  MasterPublicIP | awk -F":" '{print substr($NF, 2, length($NF)-1)}')

### Copy the credentials to cluster to avoid running config
ssh  -o "StrictHostKeyChecking no"  -i ~/.ssh/ee-cfn-virginia.pem ec2-user@${ip} "mkdir -p .aws"
scp -q -i ~/.ssh/ee-cfn-virginia.pem ~/.aws/credentials ec2-user@${ip}:.aws/
scp -q -i ~/.ssh/ee-cfn-virginia.pem ~/.aws/config ec2-user@${ip}:.aws/
scp -q -r -i ~/.ssh/ee-cfn-virginia.pem ~/MSLSP ec2-user@${ip}:/shared/code/
ssh  -o "StrictHostKeyChecking no"  -i ~/.ssh/ee-cfn-virginia.pem ec2-user@${ip} "cd /shared/code/MSLSP/AWS/;chmod 705 MSLSP_submitTiles_AWS.sh;./MSLSP_submitTiles_AWS.sh ${tileList}"



