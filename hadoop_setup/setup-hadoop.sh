#!/bin/bash

cd ~

wget http://apache.mirrors.pair.com/hadoop/common/hadoop-1.0.2/hadoop-1.0.2.tar.gz

tar xvf hadoop-1.0.2.tar.gz

mv hadoop-1.0.2 hadoop

cd ~/hadoop/conf

rm core-site.xml hdfs-site.xml hadoop-env.sh mapred-site.xml
wget http://sw.cs.wwu.edu/~bucherj/core-site.xml
wget http://sw.cs.wwu.edu/~bucherj/hadoop-env.sh
wget http://sw.cs.wwu.edu/~bucherj/hdfs-site.xml
wget http://sw.cs.wwu.edu/~bucherj/mapred-site.xml

sed -i "s/USER/$USER/g" core-site.xml

mkdir /tmp/$USER

mkdir -p ~/.ssh

pushd ~/.ssh

rm $HOME/.ssh/id_rsa_hadoop.pub
rm $HOME/.ssh/id_rsa_hadoop

echo "id_rsa_hadoop" | ssh-keygen -t rsa -P ""

cat $HOME/.ssh/id_rsa_hadoop.pub >> $HOME/.ssh/authorized_keys

popd

~/hadoop/bin/hadoop namenode -format

cd ~/hadoop/bin

wget http://sw.cs.wwu.edu/~bucherj/start-hadoop.sh
wget http://sw.cs.wwu.edu/~bucherj/stop-hadoop.sh

chmod 755 start-hadoop.sh
chmod 755 stop-hadoop.sh

mkdir ~/hadoop/test
cd ~/hadoop/test

wget http://sw.cs.wwu.edu/~bucherj/hadoop-test/test.tar.gz

tar xvf test.tar.gz
