# Deploying to Kubernetes

The important thing about deployment is that you can either deploy the perl version (which has more features) or the python version (which is probably slightly more easily understood). If you'd like the perl version, edit the container image field in the deploy to `mrbende/funce-image:01`, and for python to `ebensma/pyfunce:latest`.

## Deployment Instructions

* Launch the deployment with `kubectl apply -f kube/deploy.yaml`
* Load the data you'd like to use on Kubernetes into the PVC using `./kube/kube-load.sh <pvc> <folder-to-copy>`
* Get a shell into the pod you've created with `kubectl exec -it <pod name> -- /bin/bash`
* If you're using the python version, you'll launch into the `/app` folder, which contains some demo materials. To find what you've copied into the PVC, navigate to the `/workspaces/your-username` folder
* Run whatever FUNC-E tests you'd like, then exit the pod and save the results to your local machine using `./kube/kube-save.sh <pvc> <pvc-path-to-results>` .