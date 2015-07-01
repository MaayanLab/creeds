FROM python:2.7

# Get pip to download and install requirements:
RUN pip install pymongo Flask

# Copy the application folder inside the container
ADD . /my_application

# Expose ports
EXPOSE 5000

# Set the default directory where CMD will execute
WORKDIR /my_application

# Set the default command to execute    
# when creating a new container
CMD python app.py 5000 0.0.0.0
