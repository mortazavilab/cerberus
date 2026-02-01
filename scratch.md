Build Cerberus docker image
```bash
docker buildx create --use
docker buildx build \
    --platform linux/amd64,linux/arm64 \
    -t ghcr.io/fairliereese/cerberus:latest \
    --push .
```
