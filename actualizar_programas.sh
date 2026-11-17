#!/bin/bash

# Verifica que estás dentro de tu repositorio (ej: 1.-Logistic_Map)
if ! git rev-parse --is-inside-work-tree > /dev/null 2>&1; then
    echo "❌ Error: No estás dentro de una carpeta de Git."
    exit 1
fi

MAQUINA=$(hostname)
FECHA=$(date +"%d-%b-%Y %H:%M")

echo "🔄 Sincronizando con GitHub..."

# DESCARGA SEGURA: Si detecta un choque, frena todo inmediatamente
if ! git pull; then
    echo "⚠️ ERROR: Se detectó un conflicto. Git protegió tus archivos locales."
    echo "Por favor haz un commit de tus cambios actuales y luego un merge."
    exit 1
fi

if [[ -z $(git status -s) ]]; then
    echo "✨ Todo está actualizado. No hay nada nuevo para subir."
else
    echo "📝 Cambios detectados. Preparando respaldo..."
    git add .
    
    read -p "Mensaje (Presiona Enter para autogenerado): " MSG
    if [ -z "$MSG" ]; then
        MSG="Avances desde $MAQUINA el $FECHA"
    fi
    
    git commit -m "$MSG"
    
    if git push; then
        echo "✅ ¡Código respaldado en la nube con éxito!"
    else
        echo "❌ Error al subir a la nube."
        exit 1
    fi
fi
