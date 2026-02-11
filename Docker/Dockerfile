# Base de R específica
FROM rocker/r-ver:4.4.1

# Instalar dependencias del sistema
RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libxt-dev \
    libx11-dev \
    libglu1-mesa-dev \
    libgl1-mesa-dev \
    bash \
    python3.10 \
    python3.10-venv \
    python3-pip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Crear un directorio de trabajo
WORKDIR /usr/src/app

# Copiar el script de instalación de R y ejecutarlo
COPY install_dependencies.R .
RUN Rscript install_dependencies.R

# Copiar dependencias python
COPY install_python_dependencies.txt .
RUN python3.10 -m pip install --upgrade pip
RUN python3.10 -m pip install -r install_python_dependencies.txt

# Copiar tu pipeline al contenedor
COPY . .

# Hacer el wrapper ejecutable
RUN chmod +x Parameters_GGNC_v2.sh

# Comando por defecto para iniciar
CMD ["bash"]
