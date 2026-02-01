
Running the streamlit app

```bash
streamlit run main.py --server.maxUploadSize 2000
```

Creating Docker image
```bash
docker buildx create --use
docker buildx build \
    --platform linux/amd64,linux/arm64 \
    -t ghcr.io/fairliereese/cerberus-viewer:latest \
    --push .
```
