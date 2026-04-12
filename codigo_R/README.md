# Presentación

Este directorio contiene código R y materiales para el manual sobre análisis de datos de microarrays.

El libro tiene casi diez años, lo que, en este campo, significa que, al menos en parte, está desactualizado. 
Esto se nota especialmente en el código R en donde algunos pequeños cambios pueden hacer que algunos fragmentos de código no funcionen o se produzcan errores. 


# Errata detectada y cambios aplicados

## Instalación de paquietes de Bioconductor

- En las primeras versiones de Bioconductor la instalación se realizaba con una función `biocLite()` que se descargaba de la web de Bioconductor.
- En las versiones del año 2020 y posteriores está disponible el paquete `BiocManager` y la instalación de nuevos paquetes se realiza con `BiocManager::install`. 

## ID vs nombres de fila en las `topTable`
- Los capítulos 5 y 6 contenían referencias a un campo llamado "ID" en una `topTable` (de hecho, en el objeto `lmFit` de donde se obtiene dicha `extrajo la tabla superior`topTable`).
- Esto era así en las versiones antiguas de `limma`. Hace unos años, los desarrolladores del paquete `limma` decidieron eliminar este campo y mantener los identificadores del conjunto de sondas como nombres de fila.
- Hemos actualizado el código para eliminar referencias a esta "ID" y cambiarlas por "rownames ()"