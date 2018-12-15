FROM nvidia/cuda:9.0-base

ENV TZ=America/Edmonton
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get -y upgrade
RUN apt-get install -y cmake g++ libopenmpi-dev ssh-client

RUN apt-get install -y curl
RUN curl -L https://packages.gitlab.com/install/repositories/runner/gitlab-runner/script.deb.sh \
    | bash
COPY ./pin-gitlab-runner.pref /etc/apt/preferences.d/pin-gitlab-runner.pref
RUN apt-get install -y gitlab-runner

RUN gitlab-runner register \
    --non-interactive \
    --url "http://git.sciridge.com/" \
    --registration-token "$CI_RUNNER_REGISTRATION_TOKEN" \
    --executor "shell" \
    --description "lbflow-runner" \
    --tag-list "lbflow" \
    --run-untagged \
    --locked="false"

RUN sed -i "s/mesg n/test -t 0 \&\& mesg n/" /root/.profile

COPY ./start.sh /start.sh
CMD ["/start.sh"]
