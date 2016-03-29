#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
g++ -O3 -Wall -c -fmessage-length=0 -pthread -std=c++11 -Wl,--no-as-needed -MMD -MP -MF"${DIR}/src/threadpool11/pool.d" -MT"${DIR}/src/threadpool11/pool.d" -o "${DIR}/src/threadpool11/pool.o" "${DIR}/src/threadpool11/pool.cpp"


g++ -O3 -Wall -c -fmessage-length=0 -pthread -std=c++11 -Wl,--no-as-needed -MMD -MP -MF"${DIR}/src/threadpool11/worker.d" -MT"${DIR}/src/threadpool11/worker.d" -o "${DIR}/src/threadpool11/worker.o" "${DIR}/src/threadpool11/worker.cpp"


g++ -O3 -Wall -c -fmessage-length=0 -pthread -std=c++11 -Wl,--no-as-needed -MMD -MP -MF"${DIR}/src/main.d" -MT"${DIR}/src/main.d" -o "${DIR}/src/main.o" "${DIR}/src/main.cpp"


g++ -pthread -std=c++11 -Wl,--no-as-needed -o ./src/sigmod  ./src/threadpool11/pool.o ./src/threadpool11/worker.o  ./src/main.o   

