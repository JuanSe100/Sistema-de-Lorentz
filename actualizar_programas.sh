#!/bin/bash

# Captura el nombre de la computadora y la fecha exacta
MAQUINA=$(hostname)
FECHA=$(date +"%d-%b-%Y %H:%M")

echo "🔄 Sincronizando con GitHub..."

# 1. Siempre descarga primero para evitar conflictos
git pull

# 2. Verifica si hay archivos nuevos o modificados
if [[ -z $(git status -s) ]]; then
    echo "✨ Todo está actualizado. No hay nada nuevo para subir."
else
    echo "📝 Cambios detectados. Preparando respaldo..."
    git add .
    
    # Te da la opción de escribir un mensaje, si solo presionas Enter, usa el automático
    read -p "Mensaje (Presiona Enter para usar autogenerado): " MSG
    if [ -z "$MSG" ]; then
        MSG="Avances desde $MAQUINA el $FECHA"
    fi
    
    git commit -m "$MSG"
    git push
    echo "✅ ¡Código respaldado en la nube con éxito!"
fi
