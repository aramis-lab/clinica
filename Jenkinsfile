// Continuous Integration script for Clinica
// www.clinica.run
// Author: mauricio.diaz@inria.fr


pipeline {
  agent none 
    stages {
      stage('Build') {
        parallel {
          stage('Build in Linux') {
            agent { label 'ubuntu' }
            steps {
              echo 'Building..'
                sh 'cd clinica && ls'
            }
          }
          stage('Build Mac') {
            agent { label 'macos' }
            steps {
              echo 'Building..'
                sh 'cd clinica/ && ls'
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
