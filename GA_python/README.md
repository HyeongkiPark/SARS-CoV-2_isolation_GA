# Source code for GA_python

This code is most aimed for docker use.

First, go into "scr" repository.

Then, run the code below.

`docker image build --rm=true -t sarscov2_isolation_ga:latest .`

When the docker image is successfully made, run the code below and run&start the container.

- For Max,

`docker run --rm=true -p 8888:8888 -v $PWD:/tmp -it sarscov2_isolation_ga /bin/bash`

- For Windows,

`docker run --rm=true -p 8888:8888 -v %cd%:/tmp -it sarscov2_isolation_ga /bin/bash`

Finally, start jupyter notebook and execute the cells.

`jupyter notebook --port 8888 --ip=0.0.0.0 --allow-root --NotebookApp.token='' --notebook-dir 'tmp'`

Open the link "localhost:8888".

> **Warning**
> If another jupyter notebook is open, the link may not open successfully.
