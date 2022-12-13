pipeline {
  triggers {
      cron('0 16 * * *')
    }
  options {
    timeout(time: 6, unit: 'HOURS')
    disableConcurrentBuilds(abortPrevious: true)
  }
  agent none
  stages {
    stage('Checkout') {
      parallel {
        stage('Linux') {
          agent {
            label 'ubuntu'
          }
          environment {
            CONDA_ENV = "$WORKSPACE/env"
            CONDA_HOME = "$HOME/miniconda"
            PATH = "$HOME/.local/bin:$PATH"
            POETRY="poetry"
          }
          stages {
            stage('Build environment') {
              steps {
                sh 'echo "Agent name is ${NODE_NAME}"'
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  make env.conda
                  conda activate $CONDA_ENV
                  conda info
                  if ! command -v $POETRY &> /dev/null
                  then
                    echo "$POETRY could not be found"
                    exit
                  else
                    echo "$($POETRY --version) installed at : $(which $POETRY)"
                  fi
                '''
              }
            }
            stage('Install Clinica') {
              steps {
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate $CONDA_ENV
                  make install
                  clinica --help
                  conda list
                '''
              }
            }
            stage('PET:nonreg:notslow') {
              environment {
                PATH = "/usr/local/Modules/bin:$PATH"
                WORK_DIR = "/mnt/data/ci/working_dir_linux/PET"
                INPUT_DATA_DIR = "/mnt/data_ci"
                TMP_DIR = "/mnt/data/ci/tmp"
                }
              steps {
                catchError(buildResult: 'FAILURE', stageResult: 'UNSTABLE'){
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate $CONDA_ENV
                  source /usr/local/Modules/init/profile.sh
                  mkdir -p $WORK_DIR
                  module load clinica.all
                  cd test
                  poetry run pytest \
                    --junitxml=./test-reports/nonregression_linux_pet.xml \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_DIR \
                    --disable-warnings \
                    --timeout=0 \
                    -n 4 \
                    -m "not slow" \n
                    ./nonregression/pipelines/test_run_pipelines_pet.py
                '''
                }
              }
              post {
                always {
                  junit 'test/test-reports/*.xml'
                }
                success {
                  sh 'rm -rf ${WORK_DIR}'
                }
              }
            }
            stage('Stats:nonreg:notslow') {
              environment {
                PATH = "/usr/local/Modules/bin:$PATH"
                WORK_DIR = "/mnt/data/ci/working_dir_linux/Stats"
                INPUT_DATA_DIR = "/mnt/data_ci"
                TMP_DIR = "/mnt/data/ci/tmp"
                }
              steps {
                catchError(buildResult: 'FAILURE', stageResult: 'UNSTABLE'){
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate $CONDA_ENV
                  mkdir -p $WORK_DIR
                  source /usr/local/Modules/init/profile.sh
                  module load clinica.all
                  cd test
                  poetry run pytest \
                    --junitxml=./test-reports/nonregression_linux_stats.xml \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_DIR \
                    --disable-warnings \
                    --timeout=0 \
                    -n 4 \
                    -m "not slow" \
                    ./nonregression/pipelines/test_run_pipelines_stats.py
                '''
                }
              }
              post {
                always {
                  junit 'test/test-reports/*.xml'
                }
                success {
                  sh 'rm -rf ${WORK_DIR}'
                }
              }
            }
            stage('ML:nonreg:notslow') {
              environment {
                PATH = "/usr/local/Modules/bin:$PATH"
                WORK_DIR = "/mnt/data/ci/working_dir_linux/ML"
                INPUT_DATA_DIR = "/mnt/data_ci"
                TMP_DIR = "/mnt/data/ci/tmp"
                }
              steps {
                catchError(buildResult: 'FAILURE', stageResult: 'UNSTABLE'){
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate $CONDA_ENV
                  mkdir -p $WORK_DIR
                  source /usr/local/Modules/init/profile.sh
                  module load clinica.all
                  cd test
                  poetry run pytest \
                    --junitxml=./test-reports/nonregression_linux_ml.xml \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_DIR \
                    --disable-warnings \
                    --timeout=0 \
                    -n 4 \
                    -m "not slow" \n
                    ./nonregression/pipelines/test_run_pipelines_ml.py
                '''
                }
              }
              post {
                always {
                  junit 'test/test-reports/*.xml'
                }
                success {
                  sh 'rm -rf ${WORK_DIR}'
                }
              }
            }
            stage('Anat:nonreg:notslow') {
              environment {
                PATH = "/usr/local/Modules/bin:$PATH"
                WORK_DIR = "/mnt/data/ci/working_dir_linux/Anat"
                INPUT_DATA_DIR = "/mnt/data_ci"
                TMP_DIR = "/mnt/data/ci/tmp"
                }
              steps {
                catchError(buildResult: 'FAILURE', stageResult: 'UNSTABLE'){
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate $CONDA_ENV
                  mkdir -p $WORK_DIR
                  source /usr/local/Modules/init/profile.sh
                  module load clinica.all
                  cd test
                  poetry run pytest \
                    --junitxml=./test-reports/nonregression_linux_anat.xml \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_DIR \
                    --disable-warnings \
                    --timeout=0 \
                    -n 4 \
                    -m "not slow" \n
                    ./nonregression/pipelines/test_run_pipelines_anat.py
                '''
                }
              }
              post {
                always {
                  junit 'test/test-reports/*.xml'
                }
                success {
                  sh 'rm -rf ${WORK_DIR}'
                }
              }
            }
            stage('DWI:nonreg:notslow') {
              environment {
                PATH = "/usr/local/Modules/bin:$PATH"
                WORK_DIR = "/mnt/data/ci/working_dir_linux/DWI"
                INPUT_DATA_DIR = "/mnt/data_ci"
                TMP_DIR = "/mnt/data/ci/tmp"
                }
              steps {
                catchError(buildResult: 'FAILURE', stageResult: 'UNSTABLE'){
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate $CONDA_ENV
                  mkdir -p $WORK_DIR
                  source /usr/local/Modules/init/profile.sh
                  module load clinica.all
                  cd test
                  poetry run pytest \
                    --junitxml=./test-reports/nonregression_linux_dwi.xml \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_DIR \
                    --disable-warnings \
                    --timeout=0 \
                    -n 4 \
                    -m "not slow" \n
                    ./nonregression/pipelines/test_run_pipelines_dwi.py
                '''
                }
              }
              post {
                always {
                  junit 'test/test-reports/*.xml'
                }
                success {
                  sh 'rm -rf ${WORK_DIR}'
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
        stage('MacOS') {
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
                  conda activate $CONDA_ENV
                  conda info
                '''
              }
            }
            stage('Install Clinica') {
              steps {
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate $CONDA_ENV
                  make install
                  clinica --help
                  conda list
                '''
              }
            }
            stage("PET:nonreg:notslow") {
               environment {
                 WORK_DIR = "/Volumes/data/working_dir_mac/PET"
                 INPUT_DATA_DIR = "/Volumes/data_ci"
                 TMP_DIR = "/Volumes/data/tmp"
               }
               steps {
                 sh '''
                   source "${CONDA_HOME}/etc/profile.d/conda.sh"
                   conda activate $CONDA_ENV
                   source /usr/local/opt/modules/init/bash
                   mkdir -p $WORK_DIR
                   module load clinica.all
                   cd test
                   poetry run pytest \
                     --junitxml=./test-reports/nonregression_mac_pet.xml \
                     --verbose \
                     --working_directory=$WORK_DIR \
                     --input_data_directory=$INPUT_DATA_DIR \
                     --basetemp=$TMP_DIR \
                     --disable-warnings \
                     --timeout=0 \
                     -n 4 \
                     -m "not slow" \
                     ./nonregression/pipelines/test_run_pipelines_pet.py
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
            stage('Stats:nonreg:notslow') {
              environment {
                WORK_DIR = "/Volumes/data/working_dir_mac/Stats"
                INPUT_DATA_DIR = "/Volumes/data_ci"
                TMP_DIR = "/Volumes/data/tmp"
              }
              steps {
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate $CONDA_ENV
                  source "${BREW_PREFIX}/opt/modules/init/bash"
                  mkdir -p $WORK_DIR
                  module load clinica.all
                  cd test
                  poetry run pytest \
                    --junitxml=./test-reports/nonregression_mac_stats.xml \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_DIR \
                    --disable-warnings \
                    --timeout=0 \
                    -n 4 \
                    -m "not slow" \
                    ./nonregression/pipelines/test_run_pipelines_stats.py
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
            stage('ML:nonreg:notslow') {
              environment {
                WORK_DIR = "/Volumes/data/working_dir_mac/ML"
                INPUT_DATA_DIR = "/Volumes/data_ci"
                TMP_DIR = "/Volumes/data/tmp"
              }
              steps {
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate $CONDA_ENV
                  source "${BREW_PREFIX}/opt/modules/init/bash"
                  mkdir -p $WORK_DIR
                  module load clinica.all
                  cd test
                  poetry run pytest \
                    --junitxml=./test-reports/nonregression_mac_ml.xml \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_DIR \
                    --disable-warnings \
                    --timeout=0 \
                    -n 4 \
                    -m "not slow" \
                    ./nonregression/pipelines/test_run_pipelines_ml.py
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
            stage('Anat:nonreg:notslow') {
              environment {
                WORK_DIR = "/Volumes/data/working_dir_mac/Anat"
                INPUT_DATA_DIR = "/Volumes/data_ci"
                TMP_DIR = "/Volumes/data/tmp"
              }
              steps {
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate $CONDA_ENV
                  source "${BREW_PREFIX}/opt/modules/init/bash"
                  mkdir -p $WORK_DIR
                  module load clinica.all
                  cd test
                  poetry run pytest \
                    --junitxml=./test-reports/nonregression_mac_anat.xml \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_DIR \
                    --disable-warnings \
                    --timeout=0 \
                    -n 4 \
                    -m "not slow" \
                    ./nonregression/pipelines/test_run_pipelines_anat.py
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
          stage('DWI:nonreg:notslow') {
              environment {
                WORK_DIR = "/Volumes/data/working_dir_mac/DWI"
                INPUT_DATA_DIR = "/Volumes/data_ci"
                TMP_DIR = "/Volumes/data/tmp"
              }
              steps {
                sh '''
                  source "${CONDA_HOME}/etc/profile.d/conda.sh"
                  conda activate $CONDA_ENV
                  source "${BREW_PREFIX}/opt/modules/init/bash"
                  mkdir -p $WORK_DIR
                  module load clinica.all
                  cd test
                  poetry run pytest \
                    --junitxml=./test-reports/nonregression_mac_dwi.xml \
                    --verbose \
                    --working_directory=$WORK_DIR \
                    --input_data_directory=$INPUT_DATA_DIR \
                    --basetemp=$TMP_DIR \
                    --disable-warnings \
                    --timeout=0 \
                    -n 4 \
                    -m "not slow" \
                    ./nonregression/pipelines/test_run_pipelines_dwi.py
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
          }
          post {
            always {
              cleanWs()
            }
          }
        }
      }
    }
  }
  post {
    failure {
      mail to: 'clinica-ci@inria.fr',
        subject: "Scheduled Pipelines: ${currentBuild.fullDisplayName}",
        body: "Something is wrong with the Scheduled Pipelines ${env.BUILD_URL}"
    }
  }
}
