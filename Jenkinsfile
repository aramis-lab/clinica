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
          stage('Build in Mac') {
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
            environment {
              PATH = "$HOME/miniconda/bin:$PATH"
              }
            steps {
            echo 'Installing Clinica sources in Linux...'
            echo 'My branch name is ${BRANCH_NAME}'
            sh 'echo "My branch name is ${BRANCH_NAME}"'
            sh 'printenv'
            script {
              echo "My conda env name is clinica_env_${BRANCH_NAME}"
              }
            sh './.jenkins/scripts/launch.sh'
            sh 'printenv'
            sh 'python --version'
            }
          }
          stage('Launch in MacOS') {
            agent { label 'macos' }
            environment {
              PATH = "$HOME/miniconda3/bin:$PATH"
              }
            steps {
            echo 'Installing Clinica sources in MacOS...'
            sh 'printenv'
            sh './.jenkins/scripts/launch.sh'
            }
          }
        }
      }
      stage('Short Tests') {
        parallel {
          stage('Instantiate Linux') {
            agent { label 'ubuntu' }
            environment {
              PATH = "$HOME/miniconda/bin:/usr/local/Modules/bin:$PATH"
              CLINICA_ENV_BRANCH = "clinica_env_$BRANCH_NAME"
              WORK_DIR_LINUX = "/mnt/data/ci/working_dir_linux"
              }
            steps {
              echo 'Testing pipeline instantation...'
              sh '''
                 source activate $CLINICA_ENV_BRANCH
                 source /usr/local/Modules/init/profile.sh
                 module load clinica.all
                 cd test
                 ln -s /mnt/data/ci/data_ci_linux ./data
                 taskset -c 0-21 pytest --verbose --working_directory=$WORK_DIR_LINUX --disable-warnings --timeout=0 -n 6 -k 'test_instantiate'
                 module purge
                 source deactivate
                 '''
            }
          }
          stage('Instantiate Mac') {
            agent { label 'macos' }
            environment {
              PATH = "$HOME/miniconda3/bin:/usr/local/Cellar/modules/4.1.2/bin:$PATH"
              CLINICA_ENV_BRANCH = "clinica_env_$BRANCH_NAME"
              }
            steps {
              echo 'Testing pipeline instantation...'
              sh '''
                 source activate $CLINICA_ENV_BRANCH
                 source /usr/local/opt/modules/init/bash
                 module load clinica.all
                 cd test
                 ln -s /Volumes/data/data_ci ./data
                 pytest --verbose --disable-warnings -k 'test_instantiate'
                 module purge
                 source deactivate
                 '''
            }
          }  
          stage('Style test') {
            agent { label 'ubuntu' }
            environment {
              PATH = "$HOME/miniconda/bin:/usr/local/Modules/bin:$PATH"
              CLINICA_ENV_BRANCH = "clinica_env_$BRANCH_NAME"
              WORK_DIR_LINUX = "/mnt/data/ci/working_dir_linux"
              }
            steps {
              echo 'Testing pipeline instantation...'
              sh '''
                 source activate $CLINICA_ENV_BRANCH
                 pytest --verbose -k 'test_coding_style'
                 source deactivate
                 '''
            }
          }
        }
      }
    }
}
