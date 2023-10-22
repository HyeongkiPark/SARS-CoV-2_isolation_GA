# tips for docker use

docker image build --rm=true -t sarscov2_isolation_ga:latest .

For Max,

docker run --rm=true -p 8888:8888 -v $PWD:/tmp -it sarscov2_isolation_ga /bin/bash

For Windows,

docker run --rm=true -p 8888:8888 -v %cd%:/tmp -it sarscov2_isolation_ga /bin/bash

jupyter notebook --port 8888 --ip=0.0.0.0 --allow-root --NotebookApp.token='' --notebook-dir 'tmp'

localhost:8888

