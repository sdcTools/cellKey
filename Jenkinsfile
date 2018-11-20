pipeline  {

	agent {
		docker {
			image	'library/r-base:3.5.1-devtools-stat-1'
		}
	}

	environment {
		MAIL_TO = 'bernhard.meindl@statistik.gv.at'
		ID_JENKI_ARTI_PUBLISH = "jenkins"
		ARTIFACTORY = 'https://arti.statistik.local/artifactory'
		REPOSITORY = 'CRAN-local'
		REPOSITORY_URL = "${env.ARTIFACTORY}/${env.REPOSITORY}"
		PACKAGES = ""
		REPOSITORY_REST_URL = "${env.ARTIFACTORY}/api/cran/${env.REPOSITORY}/sources"
	}

	stages {
		stage('Publish') {
			steps {
				script {

				  // current name and version of the package
				  def name = sh (returnStdout: true, script: "cat DESCRIPTION | grep ^Package: | awk 'BEGIN {F = \":\"}; {print \$2}'").trim()
					def version = sh (returnStdout: true, script: "cat DESCRIPTION | grep ^Version: | awk 'BEGIN {F = \":\"}; {print \$2}'").trim()

				  // version of the package published on artifactory
				  def versionPublished = sh (returnStdout: true, script: "Rscript scripts-jenkins/version-artifactory.R ${env.REPOSITORY_URL}").trim()

          // Pakete nach installieren
          sh (returnStdout: true, script: "Rscript scripts-jenkins/installPackages.R ${env.PACKAGES}").trim()

          //check
          sh 'Rscript scripts-jenkins/check.R'
          //test
          sh 'Rscript scripts-jenkins/test.R'
					// package not yet in artifactory
					if (versionPublished == "NA") {

						// build package
						sh 'Rscript scripts-jenkins/build.R'

						// publish to artifactory
						withCredentials([usernamePassword(credentialsId: "${env.ID_JENKI_ARTI_PUBLISH}", passwordVariable: 'PASSWORD', usernameVariable: 'USERNAME')]) {
							sh "curl -u$USERNAME:$PASSWORD -T ${name}_${version}_R_x86_64-pc-linux-gnu.tar.gz -XPOST ${env.REPOSITORY_REST_URL}"
						}

					} else {

						echo "Version of the package: ${version}"
						echo "Version of published package: ${versionPublished}"

						// version changed
						if (version != versionPublished) {

							// build package
							sh 'Rscript scripts-jenkins/build.R'

							// publish to artifactory
							withCredentials([usernamePassword(credentialsId: "${env.ID_JENKI_ARTI_PUBLISH}", passwordVariable: 'PASSWORD', usernameVariable: 'USERNAME')]) {
								sh "curl -u$USERNAME:$PASSWORD -T ${name}_${version}_R_x86_64-pc-linux-gnu.tar.gz -XPOST ${env.REPOSITORY_REST_URL}"
							}
						} else {
							echo "Version is already published"
							currentBuild.result = 'ABORTED'
						}
					}
				}
			}
		}
	}
	post {
		success {
			emailext to: "${env.MAIL_TO}", subject: "SUCCESS: Jenkins: ${env.JOB_NAME} [${env.BUILD_NUMBER}]", body: "${env.JOB_NAME}:\n${env.BUILD_URL}"
			cleanWs()
		}
		failure {
			emailext to: "${env.MAIL_TO}", subject: "ERROR: Jenkins: ${env.JOB_NAME} [${env.BUILD_NUMBER}]", body: "${env.JOB_NAME}:\n${env.BUILD_URL}"
			cleanWs()
		}
	}
}