VERSION = '0.1'.freeze

task :default do
  system("sudo docker build --tag tnpm:#{VERSION} --rm .")
  system("sudo docker image tag tnpm:#{VERSION} tnpm:latest")
end
