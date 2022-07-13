pipeline {
  options {
    timeout(time: 2, unit: 'HOURS')
    disableConcurrentBuilds(abortPrevious: true)
  }
  agent none
  stages {
    stage('Checkout') {
      failFast true
      parallel {
        stage('Test on Linux') {
          agent {
            label 'ubuntu'
          }
          environment {
            CONDA_ENV = "$WORKSPACE/env"
            CONDA_HOME = "$HOME/miniconda"
          }
          stages {
            stage('Build environment') {
              steps {
                sh 'echo "Agent name is ${NODE_NAME}"'
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  make env.conda
                  conda activate "${CONDA_ENV}"
                  conda info
                '''
              }
            }
            stage('Install Clinica') {
              steps {
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate "${CONDA_ENV}"
                  make install
                  clinica --help
                  conda list
                '''
              }
            }
            stage('Test instantiation') {
              environment {
                PATH = "/usr/local/Modules/bin:$PATH"
                WORK_DIR = '/mnt/data/ci/working_dir_linux'
                INPUT_DATA_DIR = '/mnt/data_ci'
              }
              steps {
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate "${CONDA_ENV}"
                  source /usr/local/Modules/init/profile.sh
                  module load clinica.all
                  cd test
                  taskset -c 0-21 poetry run pytest \
                    --junitxml=./test-reports/instantiation_linux.xml \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --disable-warnings \
                    --timeout=0 \
                    -n 4 \
                    ./instantiation/
                '''
              }
              post {
                always {
                  junit 'test/test-reports/*.xml'
                }
                success {
                  sh 'rm -rf ${WORK_DIR}/*'
                }
              }
            }
            stage('Test converters') {
              environment {
                PATH = "/usr/local/Modules/bin:$PATH"
                WORK_DIR = '/mnt/data/ci/working_dir_linux'
                INPUT_DATA_DIR = '/mnt/data_ci'
                TMP_BASE = '/mnt/data/ci/tmp'
              }
              steps {
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate "${CONDA_ENV}"
                  source /usr/local/Modules/init/profile.sh
                  module load clinica.all
                  cd test
                  taskset -c 0-21 poetry run pytest \
                    --junitxml=./test-reports/run_converters_linux.xml \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_BASE \
                    --disable-warnings \
                    --timeout=0 \
                    -n 2 \
                    ./nonregression/iotools/test_run_converters.py
                  '''
              }
              post {
                always {
                  junit 'test/test-reports/*.xml'
                }
                success {
                  sh '''
                    rm -rf ${WORK_DIR}/*
                    rm -rf ${TMP_BASE}/*
                  '''
                }
              }
            }
            stage('Test iotools') {
              environment {
                PATH = "/usr/local/Modules/bin:$PATH"
                WORK_DIR = '/mnt/data/ci/working_dir_linux'
                INPUT_DATA_DIR = '/mnt/data_ci'
                TMP_BASE = '/mnt/data/ci/tmp'
              }
              steps {
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate "${CONDA_ENV}"
                  source /usr/local/Modules/init/profile.sh
                  module load clinica.all
                  cd test
                  taskset -c 0-21 poetry run pytest \
                    --junitxml=./test-reports/run_converters_linux.xml \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_BASE \
                    --disable-warnings \
                    --timeout=0 \
                    -n 2 \
                    ./nonregression/iotools/test_run_utils.py
                  '''
              }
              post {
                always {
                  junit 'test/test-reports/*.xml'
                }
                success {
                  sh '''
                    rm -rf ${WORK_DIR}/*
                    rm -rf ${TMP_BASE}/*
                  '''
                }
              }
            }
          }
          post {
            always {
              cleanWs()
            }
          }
        }
        stage('Test on macOS') {
          agent {
            label 'macos'
          }
          environment {
            CONDA_ENV = "$WORKSPACE/env"
            CONDA_HOME = "$HOME/miniconda3"
          }
          stages {
            stage('Build environment') {
              steps {
                sh 'echo "Agent name is ${NODE_NAME}"'
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  make env.conda
                  conda activate "${CONDA_ENV}"
                  conda info
                '''
              }
            }
            stage('Install Clinica') {
              steps {
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate "${CONDA_ENV}"
                  make install
                  clinica --help
                  conda list
                '''
              }
            }
            stage('Test instantiation') {
              environment {
                WORK_DIR = '/Volumes/data/working_dir_mac'
                INPUT_DATA_DIR = '/Volumes/data_ci'
              }
              steps {
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate "${CONDA_ENV}"
                  source "${BREW_PREFIX}/opt/modules/init/bash"
                  module load clinica.all
                  cd test
                  poetry run pytest \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --junitxml=./test-reports/instantation_mac.xml \
                    --disable-warnings \
                    ./instantiation/
                '''
              }
              post {
                always {
                  junit 'test/test-reports/*.xml'
                }
                success {
                  sh 'rm -rf ${WORK_DIR}/*'
                }
              }
            }
            stage('Test converters') {
              environment {
                WORK_DIR = '/Volumes/data/working_dir_mac'
                INPUT_DATA_DIR = '/Volumes/data_ci'
                TMP_BASE = '/Volumes/data/tmp'
              }
              steps {
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate "${CONDA_ENV}"
                  source "${BREW_PREFIX}/opt/modules/init/bash"
                  module load clinica.all
                  cd test
                  poetry run pytest \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_BASE \
                    --junitxml=./test-reports/run_converters_mac.xml \
                    --disable-warnings \
                    ./nonregression/iotools/test_run_converters.py
                '''
              }
              post {
                always {
                  junit 'test/test-reports/*.xml'
                }
                success {
                  sh '''
                    rm -rf ${WORK_DIR}/*
                    rm -rf ${TMP_BASE}/*
                  '''
                }
              }
            }
            stage('Test iotools') {
              environment {
                WORK_DIR = '/Volumes/data/working_dir_mac'
                INPUT_DATA_DIR = '/Volumes/data_ci'
                TMP_BASE = '/Volumes/data/tmp'
              }
              steps {
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate "${CONDA_ENV}"
                  source "${BREW_PREFIX}/opt/modules/init/bash"
                  module load clinica.all
                  cd test
                  poetry run pytest \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_BASE \
                    --junitxml=./test-reports/run_converters_mac.xml \
                    --disable-warnings \
                    ./nonregression/iotools/test_run_utils.py
                '''
              }
              post {
                always {
                  junit 'test/test-reports/*.xml'
                }
                success {
                  sh '''
                    rm -rf ${WORK_DIR}/*
                    rm -rf ${TMP_BASE}/*
                  '''
                }
              }
            }
          }
          post {
            always {
              cleanWs()
            }
          }
        }
        stage('Build and publish docs') {
          agent {
            label 'ubuntu'
          }
          environment {
            CONDA_ENV = "$WORKSPACE/env"
            CONDA_HOME = "$HOME/miniconda"
          }
          steps {
            sh '''
              source "${CONDA_HOME}/etc/profile.d/conda.sh"
              make env.doc
              conda activate "${CONDA_ENV}"
              make doc
              mv site "${CHANGE_ID:-$BRANCH_NAME}"
              scp -r "${CHANGE_ID:-$BRANCH_NAME}" aramislab:/srv/local/clinica/docs/public/
            '''
          }
          post {
            always {
              cleanWs()
            }
          }
        }
      }
    }
    stage('Deploy') {
      agent {
        label 'ubuntu'
      }
      when {
        buildingTag()
      }
      environment {
        CONDA_ENV = "$WORKSPACE/env"
        CONDA_HOME = "$HOME/miniconda"
      }
      steps {
        withCredentials(
          [
            usernamePassword(
              credentialsId: 'jenkins-token-for-pypi-clinica',
              usernameVariable: 'USERNAME',
              passwordVariable: 'PASSWORD'
            )
          ]
        ) {
          sh '''
            source "${CONDA_HOME}/etc/profile.d/conda.sh"
            make env.conda
            conda activate "${CONDA_ENV}"
            poetry publish --build -u "${USERNAME}" -p "${PASSWORD}"
          '''
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
