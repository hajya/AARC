#!/bin/bash

HADOOP_PATH=/data/hadoop/bin/bin/hadoop
HADOOP_STREAMING_PATH=/data/hadoop/bin/contrib/streaming/hadoop-*streaming*.jar 


MAPPER_FILE=./map.py
REDUCER_FILE=./reduce.py
AARC_FILE=./aarc.py
SCORING_FILE=./scoring.py


INPUT_LOCATION=/user/bucherj/aarc
OUTPUT_LOCATION=/user/bucherj/aarc-output

$HADOOP_PATH dfs -rmr $OUTPUT_LOCATION

$HADOOP_PATH jar $HADOOP_STREAMING_PATH -D mapred.max.split.size=10MB -D mapred.maps.tasks=1 -D mapred.reduce.tasks=6 -file $AARC_FILE -file $SCORING_FILE -file $MAPPER_FILE -mapper $MAPPER_FILE -file $REDUCER_FILE -reducer $REDUCER_FILE -input $INPUT_LOCATION -output $OUTPUT_LOCATION 
