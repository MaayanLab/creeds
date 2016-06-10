FROM debian:stable

RUN apt-get update

RUN apt-get -y install vim

RUN apt-get -y install python
RUN apt-get -y install python-dev
RUN apt-get -y install python-pip
RUN apt-get -y install python-setuptools

# Get pip and install numpy/scipy dependencies
RUN apt-get update && apt-get install -y build-essential gfortran libatlas-base-dev


# Copy the application folder inside the container
ADD . /creeds

RUN pip install -U pip==8.1.2

# Install required python packages
# RUN pip install -r /creeds/requirements.txt
RUN pip install numpy==1.11.0
RUN pip install scipy==0.17.0
RUN pip install pandas==0.16.0
RUN pip install requests==2.6.2
RUN pip install pymongo==3.2.2
RUN pip install joblib==0.8.4


# Install apache2 server
RUN apt-get -y install apache2 apache2-prefork-dev

RUN pip install mod_wsgi
RUN pip install -Iv Flask==0.10.1

# Expose ports
EXPOSE 80

# Set the default directory where CMD will execute
WORKDIR /creeds

# Set the default command to execute    
# when creating a new container
# CMD python run.py production
CMD bash /creeds/boot.sh
