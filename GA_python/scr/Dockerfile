FROM python:3.9

  RUN apt-get update &

  RUN pip install --upgrade pip

  COPY r_requirements.txt .

  RUN pip install -r r_requirements.txt