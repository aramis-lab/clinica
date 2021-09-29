#!/usr/bin/env groovy

pipeline {
  agent none
  stages {
    stage("Checkout") {
      failFast true
      parallel {
        stage("Test on Linux") {
          agent {
            label "ubuntu"
          }
          stages {
            stage("Build on Linux") {
              environment {
                PATH = "$HOME/miniconda/bin:$PATH"
              }
              steps {
                echo "Building environment clinica_env_${BRANCH_NAME}"
                sh '''
                  eval "$(conda shell.bash hook)"
                  conda env create --force --file environment.yml -n clinica_env_${BRANCH_NAME}
                  conda activate clinica_env_$BRANCH_NAME
                  pip install -r requirements-dev.txt
                  conda info --envs
                  '''
              }
            }
            stage("Install on Linux") {
              environment {
                PATH = "$HOME/miniconda/bin:$PATH"
              }
              steps {
                echo "Installing Clinica in clinica_env_${BRANCH_NAME}"
                sh '''
                  eval "$(conda shell.bash hook)"
                  conda activate "clinica_env_$BRANCH_NAME"
                  pip install --ignore-installed .
                  clinica --help
                  conda deactivate
                  '''
              }
            }
            stage("Test instantiation") {
              environment {
                PATH = "$HOME/miniconda/bin:/usr/local/Modules/bin:$PATH"
                CLINICA_ENV_BRANCH = "clinica_env_$BRANCH_NAME"
                WORK_DIR = "/mnt/data/ci/working_dir_linux"
                INPUT_DATA_DIR = "/mnt/data_ci"
                }
              steps {
                sh '''
                  eval "$(conda shell.bash hook)"
                  conda activate "clinica_env_$BRANCH_NAME"
                  source /usr/local/Modules/init/profile.sh
                  module load clinica.all
                  cd test
                  taskset -c 0-21 pytest \
                    --junitxml=./test-reports/instantation_linux.xml \
                    --verbose \
                    --working_directory="$WORK_DIR" \
                    --input_data_directory="$INPUT_DATA_DIR" \
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
                  junit "test/test-reports/*.xml"
                }
                success {
                  sh "rm -rf $WORK_DIR/*"
                }
              }
            }
            stage("Test converters") {
              environment {
                PATH = "$HOME/miniconda/bin:/usr/local/Modules/bin:$PATH"
                CLINICA_ENV_BRANCH = "clinica_env_$BRANCH_NAME"
                WORK_DIR = "/mnt/data/ci/working_dir_linux"
                INPUT_DATA_DIR = "/mnt/data_ci"
                TMP_BASE = "/mnt/data/ci/tmp"
              }
              steps {
                sh '''
                  eval "$(conda shell.bash hook)"
                  conda activate "clinica_env_$BRANCH_NAME"
                  source /usr/local/Modules/init/profile.sh
                  module load clinica.all
                  cd test
                  taskset -c 0-21 pytest \
                    --junitxml=./test-reports/run_converters_linux.xml \
                    --verbose \
                    --working_directory="$WORK_DIR" \
                    --input_data_directory="$INPUT_DATA_DIR" \
                    --basetemp="$TMP_BASE" \
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
                  junit "test/test-reports/*.xml"
                }
                success {
                  sh "rm -rf $WORK_DIR/*"
                  sh "rm -rf $TMP_BASE/*"
                }
              }
            }
            stage("Test iotools") {
              environment {
                PATH = "$HOME/miniconda/bin:/usr/local/Modules/bin:$PATH"
                CLINICA_ENV_BRANCH = "clinica_env_$BRANCH_NAME"
                WORK_DIR = "/mnt/data/ci/working_dir_linux"
                INPUT_DATA_DIR = "/mnt/data_ci"
                TMP_BASE = "/mnt/data/ci/tmp"
              }
              steps {
                sh '''
                  eval "$(conda shell.bash hook)"
                  conda activate "clinica_env_$BRANCH_NAME"
                  source /usr/local/Modules/init/profile.sh
                  module load clinica.all
                  cd test
                  taskset -c 0-21 pytest \
                    --junitxml=./test-reports/run_utils_linux.xml \
                    --verbose \
                    --working_directory="$WORK_DIR" \
                    --input_data_directory="$INPUT_DATA_DIR" \
                    --basetemp="$TMP_BASE" \
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
                  junit "test/test-reports/*.xml"
                }
                success {
                  sh "rm -rf $WORK_DIR/*"
                  sh "rm -rf $TMP_BASE/*"
                }
              }
            }
          }
        }
        stage("Test on macOS") {
          agent {
            label "macos"
          }
          stages {
            stage("Build on macOS") {
              environment {
                PATH = "$HOME/miniconda3/bin:$PATH"
              }
              steps {
                echo 'Building environment clinica_env_${BRANCH_NAME}'
                sh '''
                  eval "$(conda shell.bash hook)"
                  conda env create --force --file environment.yml -n clinica_env_${BRANCH_NAME}
                  conda activate clinica_env_$BRANCH_NAME
                  pip install -r requirements-dev.txt
                  conda info --envs
                 '''
              }
            }
            stage("Install on macOS") {
              environment {
                PATH = "$HOME/miniconda3/bin:$PATH"
              }
              steps {
                echo "Installing Clinica in clinica_env_${BRANCH_NAME}"
                sh '''
                  eval "$(conda shell.bash hook)"
                  conda activate "clinica_env_$BRANCH_NAME"
                  pip install --ignore-installed .
                  clinica --help
                  conda deactivate
                  '''
              }
            }
            stage("Test instantiation") {
              environment {
                PATH = "$HOME/miniconda3/bin:/usr/local/Cellar/modules/4.1.2/bin:$PATH"
                CLINICA_ENV_BRANCH = "clinica_env_$BRANCH_NAME"
                WORK_DIR = "/Volumes/data/working_dir_mac"
                INPUT_DATA_DIR = "/Volumes/data_ci"
              }
              steps {
                sh '''
                  eval "$(conda shell.bash hook)"
                  conda activate "clinica_env_$BRANCH_NAME"
                  source /usr/local/opt/modules/init/bash
                  module load clinica.all
                  cd test
                  pytest \
                    --verbose \
                    --working_directory="$WORK_DIR" \
                    --input_data_directory="$INPUT_DATA_DIR" \
                    --junitxml=./test-reports/instantation_mac.xml \
                    --disable-warnings \
                    ./instantiation/
                  module purge
                  conda deactivate
                  '''
              }
              post {
                always {
                  junit "test/test-reports/*.xml"
                }
                success {
                  sh "rm -rf $WORK_DIR/*"
                }
              }
            }
            stage("Test converters") {
              environment {
                PATH = "$HOME/miniconda3/bin:/usr/local/Cellar/modules/4.1.2/bin:$PATH"
                CLINICA_ENV_BRANCH = "clinica_env_$BRANCH_NAME"
                WORK_DIR = "/Volumes/data/working_dir_mac"
                INPUT_DATA_DIR = "/Volumes/data_ci"
                TMP_BASE = "/Volumes/data/tmp"
              }
              steps {
                sh '''
                  eval "$(conda shell.bash hook)"
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
                  junit "test/test-reports/*.xml"
                }
                success {
                  sh "rm -rf $WORK_DIR/*"
                  sh "rm -rf $TMP_BASE/*"
                }
              }
            }
            stage("Test iotools") {
              environment {
                PATH = "$HOME/miniconda3/bin:/usr/local/Cellar/modules/4.1.2/bin:$PATH"
                CLINICA_ENV_BRANCH = "clinica_env_$BRANCH_NAME"
                WORK_DIR = "/Volumes/data/working_dir_mac"
                INPUT_DATA_DIR = "/Volumes/data_ci"
                TMP_BASE = "/Volumes/data/tmp"
              }
              steps {
                sh '''
                  eval "$(conda shell.bash hook)"
                  conda activate "clinica_env_$BRANCH_NAME"
                  source /usr/local/opt/modules/init/bash
                  module load clinica.all
                  cd test
                  pytest \
                    --verbose \
                    --working_directory="$WORK_DIR" \
                    --input_data_directory="$INPUT_DATA_DIR" \
                    --basetemp="$TMP_BASE" \
                    --junitxml=./test-reports/run_utils_mac.xml \
                    --disable-warnings \
                    ./nonregression/iotools/test_run_utils.py
                  module purge
                  conda deactivate
                  '''
              }
              post {
                always {
                  junit "test/test-reports/*.xml"
                }
                success {
                  sh "rm -rf $WORK_DIR/*"
                  sh "rm -rf $TMP_BASE/*"
                }
              }
            }
          }
        }
      }
    }
  }
}