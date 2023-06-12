#!/usr/bin/env groovy

// Continuous Integration script to build mkdocs in Docker
// Author: mauricio.diaz@inria.fr

pipeline {
  agent none
  stages {
    stage('Build the docs') {
      agent {
        label 'ubuntu'
      }
      environment {
        CONDA_ENV = "$WORKSPACE/env"
        CONDA_HOME = "$HOME/miniconda3"
        PATH = "$HOME/.local/bin:$PATH"
      }
      steps {
        sh '''
          source "${CONDA_HOME}/etc/profile.d/conda.sh"
          conda create -y -p "${CONDA_ENV}" python=3.9
          conda activate "${CONDA_ENV}"
          make doc
        '''
        stash(name: 'doc_html', includes: 'site/**')
      }
      post {
        always {
          cleanWs()
        }
      }
    }
    stage('Deploy') {
      agent { label 'ubuntu' }
      steps {
        echo 'Deploying in webserver...'
        sh 'echo "Agent name: ${NODE_NAME}"'
        sh 'echo "My branch name is ${BRANCH_NAME}"'
        sh 'echo "My branch name is ${TAG_NAME}"'
        unstash(name: 'doc_html')
        sh '''
          [[ -z "${TAG_NAME}" ]] && DOCS_BRANCH="${BRANCH_NAME}" || DOCS_BRANCH="${TAG_NAME}"
          mv site "${DOCS_BRANCH}"
        '''
        sshPublisher(
          publishers: [
            sshPublisherDesc(
              configName: 'web',
              transfers: [
                sshTransfer(
                  cleanRemote: false,
                  excludes: '',
                  execCommand: '',
                  execTimeout: 120000,
                  flatten: false,
                  makeEmptyDirs: false,
                  noDefaultExcludes: false,
                  patternSeparator: '[, ]+',
                  remoteDirectory: 'clinica/docs/public/',
                  remoteDirectorySDF: false,
                  removePrefix: '',
                  sourceFiles: "${DOCS_BRANCH}/**"
                )
              ],
              usePromotionTimestamp: false,
              useWorkspaceInPromotion: false,
              verbose: false
            )
          ]
        )
        echo 'Finish uploading artifacts'
      }
    }
  }
  post {
    failure {
      mail to: 'clinica-ci@inria.fr',
        subject: "Build of Clinica's user documentation failed: ${currentBuild.fullDisplayName}",
        body: "Something went wrong during the build of Clinica's user documentation, check ${env.BUILD_URL} for details."
    }
  }
}
