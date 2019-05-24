// Continuous Integration script for Clinica
// www.clinica.run
// Author: mauricio.diaz@inria.fr

pipeline {
  agent none
    stages {
      stage('Build Env') {
        parallel {
          stage('Build in Linux') {
            agent { label 'ubuntu' }
            when { changeset "environment.yml" }
            steps {
              echo 'Building Conda environment... ${BRANCH_NAME}'
              sh 'ls'
              sh 'conda env create --force --file environment.yml -n clinica_env_${BRANCH_NAME}'
            }
          }
          stage('Build Mac') {
            agent { label 'macos' }
            when { changeset "environment.yml" }
            steps {
              echo 'Building Conda environment...' + 'env.BRANCH_NAME'
              sh 'ls'
              sh 'conda env create --force --file environment.yml -n clinica_env_${BRANCH_NAME}'
            }
          }
        }
      }
      stage('Install') {
        parallel {
          stage('Launch in Linux') {
            agent { label 'ubuntu' }
            steps {
            echo 'Installing Clinica sources in Linux...'
            echo 'My branch name is ${BRANCH_NAME}'
            sh 'echo "My branch name is ${BRANCH_NAME}"'
            sh 'printenv'
            script {
              echo "My conda env name is clinica_env_${BRANCH_NAME}"
              }
            sh './.jenkins/scripts/launch.sh'
            }
          }
          stage('Launch in MacOS') {
            agent { label 'macos' }
            steps {
            echo 'Installing Clinica sources in MacOS...'
            sh './.jenkins/scripts/launch.sh'
            }
          }
        }
      }
      stage('Test') {
        parallel {
          stage('Test Linux') {
            agent { label 'ubuntu' }
            steps {
              echo 'Testing..'
              sh 'cd clinica/ && ls'
            }
          }
          stage('Test Mac') {
            agent { label 'macos' }
            steps {
              echo 'Testing..'
              sh 'cd clinica/ && ls'
            }
          }
        }
      }
    }
}
