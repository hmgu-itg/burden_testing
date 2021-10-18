# Mummy (the wrapped MONSTER)

This is a pipeline to run genome-wide burdent tests using sequencing data. Head over to [the wiki](https://github.com/hmgu-itg/burden_testing/wiki) for detailed instructions on how to run it.



## Build
```bash
VERSION=1.5.4 # Change this appropriately 

# Build the docker image
sudo docker build \
  --build-arg BUILD_DATE="$(date -u +'%Y-%m-%dT%H:%M:%SZ')" \
  --build-arg VCS_REF="$(git rev-parse HEAD)" \
  --build-arg VERSION="$VERSION" \
  --tag burden_testing:"$VERSION" \
  --tag burden_testing:latest \
  .

sudo SINGULARITY_NOHTTPS=1 singularity build burden_testing_${VERSION} docker-daemon://burden_testing:"$VERSION"
```


sudo SINGULARITY_NOHTTPS=1 /opt/singularity371/bin/singularity build burden_testing_${VERSION} docker-daemon://burden_testing:"$VERSION"
