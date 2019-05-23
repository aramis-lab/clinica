// Continuous Integration script for Clinica
// www.clinica.run
// Author: mauricio.diaz@inria.fr

pipeline {
  agent none
    environment {
      CLINICA_ENV_BRANCH = 'clinica_env' + 'env.BRANCH_NAME'
    }
    stages {
      stage('Build Env') {
        parallel {
          stage('Build in Linux') {
            agent { label 'ubuntu' }
            when { changeset "environment.yml" }
            steps {
              echo 'Building Conda environment...' + 'env.BRANCH_NAME'
              sh 'ls'
              sh 'conda env create --force --file environment.yml -n ${env.CLINICA_ENV_BRANCH}'
            }
          }
          stage('Build Mac') {
            agent { label 'macos' }
            when { changeset "environment.yml" }
            steps {
              echo 'Building Conda environment...' + 'env.BRANCH_NAME'
              sh 'ls'
              sh 'conda env create --force --file environment.yml -n ${env.CLINICA_ENV_BRANCH}'
            }
          }
        }
      }
      stage('Install') {
        parallel {
          stage('Launch in Linux') {
            agent { label 'ubuntu' }
            steps {
            echo 'Installing Clinica in Linux...'
            echo 'My branch name is ${env.BRANCH_NAME}'
            sh 'echo "My conda env name is ${env.CLINICA_ENV_BRANCH}"'
            }
          }
          stage('Launch in MacOS') {
            agent { label 'macos' }
            steps {
            echo 'Installing Clinica in MacOS...'

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
