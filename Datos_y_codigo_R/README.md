# Presentación

Este directorio contiene código R y materiales para el manual sobre análisis de datos de microarrays.

El libro tiene casi diez años, lo que, en este campo, significa que, al menos en parte, está desactualizado.

- En una primera iteración convertiremos el documento de su formato anterior a Rmarkdown.
- Una vez hecho esto el documento se actualizará para ser consistente con el estado del arte de 2020.

## Errata detectada y cambios aplicados

### ID vs nombres de fila
- Los capítulos 5 y 6 contenían referencias a un campo llamado "ID" en una tabla superior (de hecho, en el objeto "lmFit" de donde se extrajo la tabla superior).
- Esto era así en las versiones antiguas de limma. Hace unos años, los mantenedores de limma decidieron eliminar este campo y mantener los identificadores del conjunto de sondas como nombres de fila.
- Hemos actualizado el código para eliminar referencias a esta "ID" y cambiarlas por "rownames ()"