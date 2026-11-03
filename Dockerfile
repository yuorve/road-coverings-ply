# Usamos una imagen oficial de Python basada en Debian (para usar apt-get)
FROM python:3.11-slim

# Evita que Python escriba archivos .pyc y asegura que los logs salgan directo a la terminal
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Instala las dependencias de Linux requeridas (incluyendo compiladores básicos para GDAL)
RUN apt-get update && apt-get install -y --no-install-recommends \
    mdbtools \
    libgdal-dev \
    gdal-bin \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Configura las variables de entorno necesarias para que Python encuentre GDAL
ENV CPLUS_INCLUDE_PATH=/usr/include/gdal
ENV C_INCLUDE_PATH=/usr/include/gdal

# Establece el directorio de trabajo dentro del contenedor
WORKDIR /app

# Copia los archivos de requerimientos de Python e instálalos (si tienes uno)
# Si no tienes requirements.txt, puedes comentar estas dos líneas
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copia el resto de tus archivos (.py y .html) al contenedor
COPY . .

# Comando por defecto para arrancar tu aplicación
CMD ["python", "app_main.py"]
