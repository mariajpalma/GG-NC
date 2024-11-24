# Base de R específica
FROM rocker/r-ver:4.2.2

# Instalar dependencias del sistema
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libxt-dev \
    bash \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Crear un directorio de trabajo
WORKDIR /usr/src/app

# Copiar el script de instalación de R y ejecutarlo
COPY install_dependencies.R .
RUN Rscript install_dependencies.R

# Copiar tu pipeline al contenedor
COPY . .

# Hacer el wrapper ejecutable
RUN chmod +x run_pipeline.sh

# Comando por defecto para iniciar
CMD ["bash"]
