# ❄️ Nieve España

Previsión de nieve para las estaciones de esquí españolas, con mapa de Europa tipo Windy
y tabla por cotas tipo Snow-Forecast. Se actualiza sola en GitHub, sin servidor propio.

## Cómo funciona

```
ECMWF Open Data ─┐
Open-Meteo ──────┼─► backend/build.py ─► site/data/*.json + malla ─► web (GitHub Pages)
AEMET (verifica) ┘        (GitHub Actions, 3 veces al día)
```

- **La web no llama a ningún servicio meteorológico.** Solo lee archivos ya calculados,
  así que da igual cuánta gente la use: nunca hay límites de peticiones.
- **ECMWF Open Data** (oficial, CC BY 4.0): capas de toda Europa y una previsión propia en
  cada estación a partir del perfil vertical de temperatura.
- **Open-Meteo**: ECMWF 9 km, ICON-EU y Météo-France AROME corregidos a cada cota, y
  40 escenarios de ICON para las probabilidades. Gratis solo para uso NO comercial.
- **AEMET OpenData** (oficial): cada mañana compara la previsión del día anterior con la
  precipitación medida y ajusta el peso de cada fuente.
- Si una fuente falla, se usan las demás y la web lo indica en "ver detalles".

## Puesta en marcha (sin instalar nada)

1. Crea una cuenta en github.com si no la tienes.
2. Pulsa **New repository**, nómbralo `nieve-espana`, márcalo como **Public** y créalo.
3. Pulsa **uploading an existing file** y arrastra el CONTENIDO de esta carpeta
   (`backend`, `site`, `verif`, `README.md`). Pulsa **Commit changes**.
4. La carpeta `.github` suele quedar oculta y no subirse al arrastrar. Créala a mano:
   **Add file → Create new file**, escribe como nombre
   `.github/workflows/actualizar.yml`, pega el contenido de ese archivo y guarda.
5. **Settings → Pages →** en "Source" elige **GitHub Actions**.
6. **Settings → Actions → General →** en "Workflow permissions" elige
   **Read and write permissions** y guarda.
7. (Recomendado) Clave gratuita de AEMET: pídela en
   https://opendata.aemet.es/centrodedescargas/altaUsuario (llega por email). Después:
   **Settings → Secrets and variables → Actions → New repository secret**,
   nombre `AEMET_API_KEY` y pega la clave.
8. Pestaña **Actions → Actualizar previsión → Run workflow**. Tarda unos 10 minutos.
9. La web queda en `https://TU-USUARIO.github.io/nieve-espana/`.

## Si algo falla

- Pestaña **Actions** → pulsa la ejecución en rojo → pulsa el paso en rojo y copia el texto.
- En la web, pulsa "ver detalles" (bajo el título) para ver el estado de cada fuente.

## Licencias y avisos

- ECMWF Open Data: CC BY 4.0, uso comercial permitido con atribución.
- AEMET: uso libre citando "© AEMET".
- Open-Meteo: la API gratuita es solo para uso no comercial. Para un producto comercial
  hace falta su plan de pago, o prescindir de él (la web funciona solo con ECMWF).
- Mapa base: OpenFreeMap (datos de OpenStreetMap, sin clave).
- GitHub Pages está pensado para proyectos personales, no para servicios comerciales.
