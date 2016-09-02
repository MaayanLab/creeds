FROM python:2.7

ENV CONFIG_OBJ="config.ProductionConfig"

# Get pip and install numpy/scipy dependencies
RUN apt-get update && apt-get install -y \
	build-essential \
	gfortran \
	libatlas-base-dev \
	nginx

# Update pip
RUN pip install --upgrade pip

ADD requirements.txt /my_application/requirements.txt
# Set the default directory where CMD will execute
WORKDIR /my_application
# Install required python packages
RUN pip install -r /my_application/requirements.txt

# Install gunicorn and supervisor
RUN pip install gunicorn
RUN pip install supervisor supervisor-stdout

# Expose ports
EXPOSE 80

# Copy the application folder inside the container
ADD . /my_application

# Set the default command to execute    
# when creating a new container
# CMD python run.py production
CMD ["/my_application/run.sh"]
