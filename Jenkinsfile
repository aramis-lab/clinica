#!/usr/bin/env groovy

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
            environment {
              PATH = "$HOME/miniconda/bin:$PATH"
              }
            when { 
              changeset 'pyproject.toml'
            }
            steps {
              echo 'My branch name is ${BRANCH_NAME}'
              echo 'Building Conda environment... clinica_env_${BRANCH_NAME}'
              sh '''
                 eval "$(conda shell.bash hook)"
                 conda env create --force --file environment.yml -n clinica_env_${BRANCH_NAME}
                 conda activate clinica_env_$BRANCH_NAME
                 poetry install --no-dev --no-root --extras test
                 '''
            }
          }
          /*
          stage('Build in Mac') {
            agent { label 'macos' }
            environment {
              PATH = "$HOME/miniconda3/bin:$PATH"
              }
            when {
              changeset 'requirements*'
            }
            steps {
              echo 'My branch name is ${BRANCH_NAME}'
              echo 'Building Conda environment... clinica_env_${BRANCH_NAME}'
              sh '''
                 eval "$(conda shell.bash hook)"
                 conda env create --force --file environment.yml -n clinica_env_${BRANCH_NAME}
                 conda activate clinica_env_$BRANCH_NAME
                 poetry install --no-dev --no-root --extras test
                 '''
            }
          }
          */
        }
      }
      stage('Install') {
        parallel {
          stage('Launch') {
            agent { label 'ubuntu' }
            environment {
              PATH = "$HOME/miniconda/bin:$PATH"
              }
            steps {
            echo 'Installing Clinica sources in Linux...'
            echo 'My branch name is ${BRANCH_NAME}'
            echo "My conda env name is clinica_env_${BRANCH_NAME}"
            sh 'printenv'
            sh 'echo "Agent name: ${NODE_NAME}"'
            sh '''
               set +x
               source ./.jenkins/scripts/find_env.sh
               conda info --envs
               eval "$(conda shell.bash hook)"
               conda activate clinica_env_$BRANCH_NAME
               echo "Install clinica using poetry..."
               poetry install --no-dev
               eval "$(register-python-argcomplete clinica)"
               # Show clinica help message
               echo "Display clinica help message"
               clinica --help
               conda deactivate
               '''
            }
          }
          /*
          stage('Launch in MacOS') {
            agent { label 'macos' }
            environment {
              PATH = "$HOME/miniconda3/bin:$PATH"
              }
            steps {
            echo 'Installing Clinica sources in MacOS...'
            sh 'echo "Agent name: ${NODE_NAME}"'
            echo 'My branch name is ${BRANCH_NAME}'
            echo "My conda env name is clinica_env_${BRANCH_NAME}"
            sh '''
               set +x
               eval "$(conda shell.bash hook)"
               echo $CONDA_PREFIX
               source ./.jenkins/scripts/find_env.sh
               conda activate clinica_env_$BRANCH_NAME
               echo "Install clinica using poetry..."
               poetry install --no-dev
               eval "$(register-python-argcomplete clinica)"
               # Show clinica help message
               echo "Display clinica help message"
               clinica --help
               conda deactivate
            '''
            }
          }
          */
        }
      }
      stage('Instantiate Tests') {
        parallel {
          stage('Instantiate Linux') {
            agent { label 'ubuntu' }
            environment {
              PATH = "$HOME/miniconda/bin:/usr/local/Modules/bin:$PATH"
              CLINICA_ENV_BRANCH = "clinica_env_$BRANCH_NAME"
              WORK_DIR = "/mnt/data/ci/working_dir_linux"
              INPUT_DATA_DIR = "/mnt/data_ci"
              }
            steps {
              echo 'Testing pipeline instantation...'
              sh 'echo "Agent name: ${NODE_NAME}"'
              sh '''
                 set +x
                 eval "$(conda shell.bash hook)"
                 source ./.jenkins/scripts/find_env.sh
                 conda activate clinica_env_$BRANCH_NAME
                 source /usr/local/Modules/init/profile.sh
                 module load clinica.all
                 cd test
                 taskset -c 0-21 pytest \
                    --junitxml=./test-reports/instantation_linux.xml \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --disable-warnings \
                    --timeout=0 \
                    -n 4 \
                   ./instantiation/
                 module purge
                 conda deactivate
                 '''
            }
            post {
              always {
                junit 'test/test-reports/*.xml'
              }
              success {
                sh '''
                   rm -rf WORK_DIR/*
                   ''' 
              }
            }  
          }
          /*
          stage('Instantiate Mac') {
            agent { label 'macos' }
            environment {
              PATH = "$HOME/miniconda3/bin:/usr/local/Cellar/modules/4.1.2/bin:$PATH"
              CLINICA_ENV_BRANCH = "clinica_env_$BRANCH_NAME"
              WORK_DIR = "/Volumes/data/working_dir_mac"
              INPUT_DATA_DIR = "/Volumes/data_ci"
              }
            steps {
              echo 'Testing pipeline instantation...'
              sh 'echo "Agent name: ${NODE_NAME}"'
              sh '''
                 set +x
                 eval "$(conda shell.bash hook)"
                 source ./.jenkins/scripts/find_env.sh
                 conda activate clinica_env_$BRANCH_NAME
                 source /usr/local/opt/modules/init/bash
                 module load clinica.all
                 cd test
                 pytest \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --junitxml=./test-reports/instantation_mac.xml \
                    --disable-warnings \
                    ./instantiation/
                 module purge
                 conda deactivate
                 '''
            }
            post {
              always {
                junit 'test/test-reports/*.xml'
              }
              success {
                sh '''
                   rm -rf WORK_DIR/*
                   ''' 
              }
            }  
          }
          */
        }
      }
      stage('Non-regression Tests') {
        parallel {
          stage('Converters Linux') {
            agent { label 'ubuntu' }
            environment {
              PATH = "$HOME/miniconda/bin:/usr/local/Modules/bin:$PATH"
              CLINICA_ENV_BRANCH = "clinica_env_$BRANCH_NAME"
              WORK_DIR = "/mnt/data/ci/working_dir_linux"
              INPUT_DATA_DIR = "/mnt/data_ci"
              TMP_BASE = "/mnt/data/ci/tmp"
              }
            steps {
              echo 'Testing pipeline instantation...'
              sh 'echo "Agent name: ${NODE_NAME}"'
              sh '''
                 set +x
                 eval "$(conda shell.bash hook)"
                 source ./.jenkins/scripts/find_env.sh
                 conda activate clinica_env_$BRANCH_NAME
                 source /usr/local/Modules/init/profile.sh
                 module load clinica.all
                 cd test
                 taskset -c 0-21 pytest \
                    --junitxml=./test-reports/run_converters_linux.xml \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_BASE \
                    --disable-warnings \
                    --timeout=0 \
                    -n 2 \
                    ./nonregression/iotools/test_run_converters.py
                 module purge
                 conda deactivate
                 '''
            }
            post {
              always {
                junit 'test/test-reports/*.xml'
              }
              success {
                sh '''
                   rm -rf WORK_DIR/*
                   rm -rf TMP_BASE/*
                   ''' 
              }
            }  
          }
          stage('Iotools Linux') {
            agent { label 'ubuntu' }
            environment {
              PATH = "$HOME/miniconda/bin:/usr/local/Modules/bin:$PATH"
              CLINICA_ENV_BRANCH = "clinica_env_$BRANCH_NAME"
              WORK_DIR = "/mnt/data/ci/working_dir_linux"
              INPUT_DATA_DIR = "/mnt/data_ci"
              TMP_BASE = "/mnt/data/ci/tmp"
              }
            steps {
              echo 'Testing pipeline instantation...'
              sh 'echo "Agent name: ${NODE_NAME}"'
              sh '''
                 set +x
                 eval "$(conda shell.bash hook)"
                 source ./.jenkins/scripts/find_env.sh
                 conda activate clinica_env_$BRANCH_NAME
                 source /usr/local/Modules/init/profile.sh
                 module load clinica.all
                 cd test
                 taskset -c 0-21 pytest \
                    --junitxml=./test-reports/run_utils_linux.xml \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_BASE \
                    --disable-warnings \
                    --timeout=0 \
                    -n 2 \
                    ./nonregression/iotools/test_run_utils.py
                 module purge
                 conda deactivate
                 '''
            }
            post {
              always {
                junit 'test/test-reports/*.xml'
              }
              success {
                sh '''
                   rm -rf WORK_DIR/*
                   rm -rf TMP_BASE/*
                   ''' 
              }
            }  
          }
          stage('Converters Mac') {
            agent { label 'macos' }
            environment {
              PATH = "$HOME/miniconda3/bin:/usr/local/Cellar/modules/4.1.2/bin:$PATH"
              CLINICA_ENV_BRANCH = "clinica_env_$BRANCH_NAME"
              WORK_DIR = "/Volumes/data/working_dir_mac"
              INPUT_DATA_DIR = "/Volumes/data_ci"
              TMP_BASE = "/Volumes/data/tmp"
              }
            steps {
              echo 'Testing pipeline instantation...'
              sh 'echo "Agent name: ${NODE_NAME}"'
              sh '''
                 set +x
                 eval "$(conda shell.bash hook)"
                 source ./.jenkins/scripts/find_env.sh
                 conda activate clinica_env_$BRANCH_NAME
                 source /usr/local/opt/modules/init/bash
                 module load clinica.all
                 cd test
                 pytest \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_BASE \
                    --junitxml=./test-reports/run_converters_mac.xml \
                    --disable-warnings \
                    ./nonregression/iotools/test_run_converters.py
                 module purge
                 conda deactivate
                 '''
            }
            post {
              always {
                junit 'test/test-reports/*.xml'
              }
              success {
                sh '''
                   rm -rf WORK_DIR/*
                   rm -rf TMP_BASE/*
                   ''' 
              }
            }  
          }
          stage('Iotools Mac') {
            agent { label 'macos' }
            environment {
              PATH = "$HOME/miniconda3/bin:/usr/local/Cellar/modules/4.1.2/bin:$PATH"
              CLINICA_ENV_BRANCH = "clinica_env_$BRANCH_NAME"
              WORK_DIR = "/Volumes/data/working_dir_mac"
              INPUT_DATA_DIR = "/Volumes/data_ci"
              TMP_BASE = "/Volumes/data/tmp"
              }
            steps {
              echo 'Testing pipeline instantation...'
              sh 'echo "Agent name: ${NODE_NAME}"'
              sh '''
                 set +x
                 eval "$(conda shell.bash hook)"
                 source ./.jenkins/scripts/find_env.sh
                 conda activate clinica_env_$BRANCH_NAME
                 source /usr/local/opt/modules/init/bash
                 module load clinica.all
                 cd test
                 pytest \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_BASE \
                    --junitxml=./test-reports/run_utils_mac.xml \
                    --disable-warnings \
                    ./nonregression/iotools/test_run_utils.py
                 module purge
                 conda deactivate
                 '''
            }
            post {
              always {
                junit 'test/test-reports/*.xml'
              }
              success {
                sh '''
                   rm -rf WORK_DIR/*
                   rm -rf TMP_BASE/*
                   ''' 
              }
            }  
          }
        }
      }
      stage('Build and publish doc') {
        agent { label 'ubuntu && short' }
        environment {
          PATH = "$HOME/miniconda/bin:$PATH"
          }
        steps {
          echo "Build and publish Clinica's documentation"
          sh 'echo "Agent name: ${NODE_NAME}"'
          sh '''
          eval "$(conda shell.bash hook)"
          conda create --yes --name build_doc python=3.8 poetry
          conda activate build_doc
          poetry install --no-dev --no-root --extras docs
          ./.jenkins/scripts/publish.sh ${BRANCH_NAME}
          '''
        }
        post {
          success {
            sh 'scp -r ${BRANCH_NAME} aramislab:~/clinica/docs/public/'
          }
        }
      }
      stage('Deploy') {
        agent { label 'ubuntu' }
        when { buildingTag() }
        environment {
          PATH = "$HOME/miniconda/bin:$PATH"
        }
        steps {
          echo 'Create ClinicaDL package and upload to Pypi...'
          sh 'echo "Agent name: ${NODE_NAME}"'
          sh '''#!/usr/bin/env bash
             set +x
             eval "$(conda shell.bash hook)"
             source ./.jenkins/scripts/find_env.sh
             conda activate clinica_env_$BRANCH_NAME
             poetry install --no-dev
             clinica --help
             cd $WORKSPACE/.jenkins/scripts
             ./generate_wheels.sh
             conda deactivate
             '''
          withCredentials([usernamePassword(credentialsId: 'jenkins-pass-for-pypi-aramis', usernameVariable: 'USERNAME', passwordVariable: 'PASSWORD')]) {
             sh '''#!/usr/bin/env bash
                cd $WORKSPACE
                twine upload \
                  -u ${USERNAME} \
                  -p ${PASSWORD} ./dist/*
                '''
             }
        }
        post {
          success {
            mattermostSend( 
              color: "#00B300",
              message: "Clinica version ${env.TAG_NAME} is now available!:  ${env.JOB_NAME} #${env.BUILD_NUMBER} (<${env.BUILD_URL}|Link to build>)"
            )
          }
        }
      }
    }
    post {
      failure {
        mail to: 'clinica-ci@inria.fr',
          subject: "Failed Pipeline: ${currentBuild.fullDisplayName}",
          body: "Something is wrong with ${env.BUILD_URL}"
      }
    }
  }
