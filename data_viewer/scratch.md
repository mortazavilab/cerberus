
Running the streamlit app

```bash
streamlit run main.py --server.maxUploadSize 2000
```

Creating Docker image (VIEWER)
```bash
cd  ~/mortazavi_lab/bin/cerberus/data_viewer
docker buildx create --use
docker buildx build \
    --platform linux/amd64,linux/arm64 \
    -t ghcr.io/fairliereese/cerberus-viewer:latest \
    --push .
```

Creating Docker image (CERBERUS ONLY)
```bash
cd  ~/mortazavi_lab/bin/cerberus/
docker buildx create --use
docker buildx build \
    --platform linux/amd64,linux/arm64 \
    -t ghcr.io/fairliereese/cerberus:latest \
    --push .
```
