# Clinica Docker Image

## Introduction

This repository contains the Clinica software Docker image files. This image embed Clinica and _all_ its dependancies in an Ubuntu.

## Getting started

### Prerequisites

- **Docker**: Obviously Docker is needed to use this Docker image. Install Docker software as you would install any general public software by [downloading](https://www.docker.com/products/overview) and executing the right installer for you operating system.
- **GitLab account**: Clinica's repository is still under restricted access, so you need a valid GitLab login and password*-** to clone it.
- **MATLAB(r) license file**: To activate a MATLAB(r) distribution on the Docker image, you will need your license file. You can download it from the [MathWorks(r) License Center](https://www.mathworks.com/licensecenter/licenses), after log in click on your license > Install and Activate > Update License File > Activate a Computer. **Important**: You must specify `root` as a Computer Login Name when you generate the license file.
- **FreeSurfer license file**: To activate your FreeSurfer software, again, a license is needed. You just need to fill in [this form](https://surfer.nmr.mgh.harvard.edu/registration.html); a mail containing the license information will be sent to you.
<!-- ARAMIS LAB users have to be connected through Wi-fi; FreeSurfer download are somehow blocked by the wired network? -->

### Installation

The first thing to do is to clone this repository:

```bash
git clone git@gitlab.icm-institute.org:aramis/clinica-docker.git
```

<!-- See details on [this link](https://docs.docker.com/engine/reference/commandline/build/#/git-repositories) about how to avoid `clone` + `build`. -->

Copy your previously downloaded MATLAB(r) license file in the `clinica-docker/` directory. Rename it to `license.lic` if it is needed:

```bash
cp <PATH_TO_YOUR_MATLAB_LICENSE_FILE> clinica-docker/license.lic
```

Paste the FreeSurfer license information inside the empty `clinica-docker/freesurfer_license.txt` file. And finally build the Clinica Docker image*:

```bash
docker build --build-arg GITLAB_PSWD='<YOUR_GITLAB_PASSWORD>' --build-arg GITLAB_USER='<YOUR_GITLAB_USERNAME>' --tag aramislab/clinica clinica-docker/
```

### Usage

To run Clinica, just type the following:

```bash
docker run aramislab/clinica
```

You then use this command as if it was `clinica` only. For instance:

```bash
docker run aramislab/clinica run mypipeline <INPUT_DIR> <OUTPUT_DIR>
```

To access data stored on your computer from a Docker image, you must share it using the `-v` option as follow:

```bash
docker run -v <HOST_DIR>:<CONTAINER_DIR> aramislab/clinica run mypipeline <INPUT_DIR> <OUTPUT_DIR>
```


## Docker basics

### Cheat sheet

Command                            | Purpose
---------------------------------- | ------------------------------------
`docker -h`                        | Display Docker help
`docker images`                    | Display the downloaded images
`docker ps [-a]`                   | Display the running [or all] containers
`docker rmi <IMAGE_NAME_OR_ID>`    | Remove an image
`docker rm <CONTAINER_NAME_OR_ID>` | Remove a container

One will often want to remove all the containers that previously ran. You can do so by using this custom bash function to add in your `~/.bashrc` file:

```bash
function docker-rm-all {
    docker stop $(docker ps -a -q) && docker rm $(docker ps -a -q)
}
```

And eventually remove all docker images:

```bash
function docker-rmi-all {
    docker rmi -f $(docker images -q)
}
```

> **Warning**: Please be aware that the latter command will remove _all_ images that you may have downloaded.

### Customizing your Clinica Docker image

By default, each time you run an image Docker will create what is called a _container_. If you want to modify your Clinica Docker image, you can do it by running it interactively:

```bash
docker start <CONTAINER_NAME_OR_ID>
docker exec -it <CONTAINER_NAME_OR_ID> bash
```

This will start an interactive bash shell from your Docker container. You can do any modification on it and then use the following command to "save" your modifications:

```bash
docker commit <CONTAINER_NAME_OR_ID> myclinica
```

See [here](https://docs.docker.com/engine/getstarted/step_six/) for more details on how to manage Docker containers and images.

> *: Your GitLab account password has to be written in clear directly in the terminal when entering the `docker build` command. Be careful not to let anyone see your command history. For more information on how to delete your terminal command history see [this page](http://unix.stackexchange.com/questions/203290/how-do-i-clear-the-terminal-history).

> **: Your GitLab username (and you password) must not contain any special character and should not be your ICM mail address. The username can be found and chosen on [this page](https://gitlab.icm-institute.org/profile/account).
