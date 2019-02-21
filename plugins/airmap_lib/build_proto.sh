#!/bin/sh

protoc --proto_path=../proto --python_out=. ../proto/messages.proto
