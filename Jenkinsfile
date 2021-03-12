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
              changeset 'requirements*'
            }
            steps {
              echo 'My branch name is ${BRANCH_NAME}'
              echo 'Building Conda environment... clinica_env_${BRANCH_NAME}'
              sh '''
                 eval "$(conda shell.bash hook)"
                 conda env create --force --file environment.yml -n clinica_env_${BRANCH_NAME}
                 conda activate clinica_env_$BRANCH_NAME
                 pip install -r requirements-dev.txt
                 '''
            }
          }
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
                 pip install -r requirements-dev.txt
                 '''
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
            echo "My conda env name is clinica_env_${BRANCH_NAME}"
            sh 'printenv'
            sh 'echo "Agent name: ${NODE_NAME}"'
            sh '''
               set +x
               source ./.jenkins/scripts/find_env.sh
               conda info --envs
               eval "$(conda shell.bash hook)"
               conda activate clinica_env_$BRANCH_NAME
               echo "Install clinica using pip..."
               pip install --ignore-installed .
               eval "$(register-python-argcomplete clinica)"
               # Show clinica help message
               echo "Display clinica help message"
               clinica --help
               conda deactivate
               '''
            }
          }
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
               echo "Install clinica using pip..."
               pip install --ignore-installed .
               eval "$(register-python-argcomplete clinica)"
               # Show clinica help message
               echo "Display clinica help message"
               clinica --help
               conda deactivate
            '''
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
              sh 'echo "Agent name: ${NODE_NAME}"'
              sh '''
                 set +x
                 eval "$(conda shell.bash hook)"
                 source ./.jenkins/scripts/find_env.sh
                 conda activate clinica_env_$BRANCH_NAME
                 source /usr/local/Modules/init/profile.sh
                 module load clinica.all
                 cd test
                 ln -s /mnt/data/ci/data_ci_linux ./data
                 taskset -c 0-21 pytest \
                    --junitxml=./test-reports/instantation_linux.xml \
                    --verbose \
                    --working_directory=$WORK_DIR_LINUX \
                    --disable-warnings \
                    --timeout=0 \
                    -n 6 \
                    -k 'test_instantiate'
                 module purge
                 conda deactivate
                 '''
            }
            post {
              always {
                junit 'test/test-reports/*.xml'
              }
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
              sh 'echo "Agent name: ${NODE_NAME}"'
              sh '''
                 set +x
                 eval "$(conda shell.bash hook)"
                 source ./.jenkins/scripts/find_env.sh
                 conda activate clinica_env_$BRANCH_NAME
                 source /usr/local/opt/modules/init/bash
                 module load clinica.all
                 cd test
                 ln -s /Volumes/data/data_ci ./data
                 pytest \
                    --verbose \
                    --junitxml=./test-reports/instantation_mac.xml \
                    --disable-warnings \
                    -k 'test_instantiate'
                 module purge
                 conda deactivate
                 '''
            }
            post {
              always {
                junit 'test/test-reports/*.xml'
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
          conda create --name build_doc python=3.8
          conda activate build_doc
          pip install -r docs/requirements.txt
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
             pip install -e .
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
