#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/ioctl.h>
#include <poll.h>
#include <array>
#include <random>
#include <list>
#include <chrono>
#include <omp.h>
#include <limits>
#include <algorithm>
#include <map>
#include <thread>
#include <csignal>
using namespace std;
using namespace std::chrono;

constexpr int N{2};//Number of players
constexpr bool Debug_AI{false},Timeout{false};
constexpr int PIPE_READ{0},PIPE_WRITE{1};
constexpr double FirstTurnTime{1*(Timeout?1:10)},TimeLimit{0.1*(Timeout?1:10)};

bool stop{false};//Global flag to stop all arena threads when SIGTERM is received

enum location{SAMPLES=0,DIAGNOSIS=1,MOLECULES=2,LABORATORY=3,START=4};
enum action_type{GOTO=0,CONNECT=1};

const vector<string> locationToStr{"SAMPLES","DIAGNOSIS","MOLECULES","LABORATORY","START_POS"};
const map<string,location> StrToLocation{{"SAMPLES",SAMPLES},{"DIAGNOSIS",DIAGNOSIS},{"MOLECULES",MOLECULES},{"LABORATORY",LABORATORY}};
const array<location,4> intToLocation{SAMPLES,DIAGNOSIS,MOLECULES,LABORATORY};
const vector<string> typeToStr{"A","B","C","D","E"};
const map<string,int> TypeToInt{{"A",0},{"B",1},{"C",2},{"D",3},{"E",4}};
constexpr array<array<int,4>,5> Distances{array<int,4>{0,3,3,3},array<int,4>{3,0,3,4},array<int,4>{3,3,0,3},array<int,4>{3,4,3,0},array<int,4>{2,2,2,2}};//Distances from everywhere to everywhere

struct sample_template{
    array<int,5> Cost;
    int score,exp;
};

struct sample:public sample_template{
    bool diagnosed;
    int id,rank,owner;
    inline void operator=(const sample_template &a)noexcept{
        Cost=a.Cost;
        score=a.score;
        exp=a.exp;
    }
    inline sample(const sample_template &a)noexcept{
        *this=a;
    }
};

const array<vector<sample_template>,3> SampleList{
    vector<sample_template>{
        {{ 0, 3, 0, 0, 0 }, 1, 0},
        {{ 0, 0, 0, 2, 1 }, 1, 0},
        {{ 0, 1, 1, 1, 1 }, 1, 0},
        {{ 0, 2, 0, 0, 2 }, 1, 0},
        {{ 0, 0, 4, 0, 0 }, 10, 0},
        {{ 0, 1, 2, 1, 1 }, 1, 0},
        {{ 0, 2, 2, 0, 1 }, 1, 0},
        {{ 3, 1, 0, 0, 1 }, 1, 0},
        {{ 1, 0, 0, 0, 2 }, 1, 1},
        {{ 0, 0, 0, 0, 3 }, 1, 1},
        {{ 1, 0, 1, 1, 1 }, 1, 1},
        {{ 0, 0, 2, 0, 2 }, 1, 1},
        {{ 0, 0, 0, 4, 0 }, 10, 1},
        {{ 1, 0, 1, 2, 1 }, 1, 1},
        {{ 1, 0, 2, 2, 0 }, 1, 1},
        {{ 0, 1, 3, 1, 0 }, 1, 1},
        {{ 2, 1, 0, 0, 0 }, 1, 2},
        {{ 0, 0, 0, 3, 0 }, 1, 2},
        {{ 1, 1, 0, 1, 1 }, 1, 2},
        {{ 0, 2, 0, 2, 0 }, 1, 2},
        {{ 0, 0, 0, 0, 4 }, 10, 2},
        {{ 1, 1, 0, 1, 2 }, 1, 2},
        {{ 0, 1, 0, 2, 2 }, 1, 2},
        {{ 1, 3, 1, 0, 0 }, 1, 2},
        {{ 0, 2, 1, 0, 0 }, 1, 3},
        {{ 3, 0, 0, 0, 0 }, 1, 3},
        {{ 1, 1, 1, 0, 1 }, 1, 3},
        {{ 2, 0, 0, 2, 0 }, 1, 3},
        {{ 4, 0, 0, 0, 0 }, 10, 3},
        {{ 2, 1, 1, 0, 1 }, 1, 3},
        {{ 2, 0, 1, 0, 2 }, 1, 3},
        {{ 1, 0, 0, 1, 3 }, 1, 3},
        {{ 0, 0, 2, 1, 0 }, 1, 4},
        {{ 0, 0, 3, 0, 0 }, 1, 4},
        {{ 1, 1, 1, 1, 0 }, 1, 4},
        {{ 2, 0, 2, 0, 0 }, 1, 4},
        {{ 0, 4, 0, 0, 0 }, 10, 4},
        {{ 1, 2, 1, 1, 0 }, 1, 4},
        {{ 2, 2, 0, 1, 0 }, 1, 4},
        {{ 0, 0, 1, 3, 1 }, 1, 4}
    },
    vector<sample_template>{
        {{ 0, 0, 0, 5, 0 }, 20, 0},
        {{ 6, 0, 0, 0, 0 }, 30, 0},
        {{ 0, 0, 3, 2, 2 }, 10, 0},
        {{ 0, 0, 1, 4, 2 }, 20, 0},
        {{ 2, 3, 0, 3, 0 }, 10, 0},
        {{ 0, 0, 0, 5, 3 }, 20, 0},
        {{ 0, 5, 0, 0, 0 }, 20, 1},
        {{ 0, 6, 0, 0, 0 }, 30, 1},
        {{ 0, 2, 2, 3, 0 }, 10, 1},
        {{ 2, 0, 0, 1, 4 }, 20, 1},
        {{ 0, 2, 3, 0, 3 }, 20, 1},
        {{ 5, 3, 0, 0, 0 }, 20, 1},
        {{ 0, 0, 5, 0, 0 }, 20, 2},
        {{ 0, 0, 6, 0, 0 }, 30, 2},
        {{ 2, 3, 0, 0, 2 }, 10, 2},
        {{ 3, 0, 2, 3, 0 }, 10, 2},
        {{ 4, 2, 0, 0, 1 }, 20, 2},
        {{ 0, 5, 3, 0, 0 }, 20, 2},
        {{ 5, 0, 0, 0, 0 }, 20, 3},
        {{ 0, 0, 0, 6, 0 }, 30, 3},
        {{ 2, 0, 0, 2, 3 }, 10, 3},
        {{ 1, 4, 2, 0, 0 }, 20, 3},
        {{ 0, 3, 0, 2, 3 }, 10, 3},
        {{ 3, 0, 0, 0, 5 }, 20, 3},
        {{ 0, 0, 0, 0, 5 }, 20, 4},
        {{ 0, 0, 0, 0, 6 }, 30, 4},
        {{ 3, 2, 2, 0, 0 }, 10, 4},
        {{ 0, 1, 4, 2, 0 }, 20, 4},
        {{ 3, 0, 3, 0, 2 }, 10, 4},
        {{ 0, 0, 5, 3, 0 }, 20, 4}
    },
    vector<sample_template>{
        {{ 0, 0, 0, 0, 7 }, 40, 0},
        {{ 3, 0, 0, 0, 7 }, 50, 0},
        {{ 3, 0, 0, 3, 6 }, 40, 0},
        {{ 0, 3, 3, 5, 3 }, 30, 0},
        {{ 7, 0, 0, 0, 0 }, 40, 1},
        {{ 7, 3, 0, 0, 0 }, 50, 1},
        {{ 6, 3, 0, 0, 3 }, 40, 1},
        {{ 3, 0, 3, 3, 5 }, 30, 1},
        {{ 0, 7, 0, 0, 0 }, 40, 2},
        {{ 0, 7, 3, 0, 0 }, 50, 2},
        {{ 3, 6, 3, 0, 0 }, 40, 2},
        {{ 5, 3, 0, 3, 3 }, 30, 2},
        {{ 0, 0, 7, 0, 0 }, 40, 3},
        {{ 0, 0, 7, 3, 0 }, 50, 3},
        {{ 0, 3, 6, 3, 0 }, 40, 3},
        {{ 3, 5, 3, 0, 3 }, 30, 3},
        {{ 0, 0, 0, 7, 0 }, 40, 4},
        {{ 0, 0, 0, 7, 3 }, 50, 4},
        {{ 0, 0, 3, 6, 3 }, 40, 4},
        {{ 3, 3, 5, 3, 0 }, 30, 4}
    }
};

struct player{
    location r;
    int eta,score;
    array<int,5> Mol,Exp;
    vector<sample> Samp;
};

struct project{
    array<int,5> Target;
};

const vector<project> ProjectList{
    { 3, 3, 0, 0, 3 },
    { 0, 3, 3, 3, 0 },
    { 3, 0, 0, 3, 3 },
    { 0, 0, 4, 4, 0 },
    { 0, 4, 4, 0, 0 },
    { 0, 0, 0, 4, 4 },
    { 4, 0, 0, 0, 4 },
    { 3, 3, 3, 0, 0 },
    { 0, 0, 3, 3, 3 },
    { 4, 4, 0, 0, 0 }
};

struct state{
    int SampleCount;
    array<player,2> P;
    array<int,5> Avail;
    vector<project> Proj;
    vector<sample> Samp;
    array<list<sample_template>,3> SamplePool;
};

struct action{
    action_type type;
    int id;
};

inline string EmptyPipe(const int fd){
    int nbytes;
    if(ioctl(fd,FIONREAD,&nbytes)<0){
        throw(4);
    }
    string out;
    out.resize(nbytes);
    if(read(fd,&out[0],nbytes)<0){
        throw(4);
    }
    return out;
}

struct AI{
    int id,pid,outPipe,errPipe,inPipe,turnOfDeath;
    string name;
    inline void stop(const int turn=-1){
        if(alive()){
            kill(pid,SIGTERM);
            int status;
            waitpid(pid,&status,0);//It is necessary to read the exit code for the process to stop
            if(!WIFEXITED(status)){//If not exited normally try to "kill -9" the process
                kill(pid,SIGKILL);
            }
            turnOfDeath=turn;
        }
    }
    inline bool alive()const{
        return kill(pid,0)!=-1;//Check if process is still running
    }
    inline void Feed_Inputs(const string &inputs){
        if(write(inPipe,&inputs[0],inputs.size())!=inputs.size()){
            throw(5);
        }
    }
    inline ~AI(){
        close(errPipe);
        close(outPipe);
        close(inPipe);
        stop();
    }
};

void StartProcess(AI &Bot){
    int StdinPipe[2];
    int StdoutPipe[2];
    int StderrPipe[2];
    if(pipe(StdinPipe)<0){
        perror("allocating pipe for child input redirect");
    }
    if(pipe(StdoutPipe)<0){
        close(StdinPipe[PIPE_READ]);
        close(StdinPipe[PIPE_WRITE]);
        perror("allocating pipe for child output redirect");
    }
    if(pipe(StderrPipe)<0){
        close(StderrPipe[PIPE_READ]);
        close(StderrPipe[PIPE_WRITE]);
        perror("allocating pipe for child stderr redirect failed");
    }
    int nchild{fork()};
    if(nchild==0){//Child process
        if(dup2(StdinPipe[PIPE_READ],STDIN_FILENO)==-1){// redirect stdin
            perror("redirecting stdin");
            return;
        }
        if(dup2(StdoutPipe[PIPE_WRITE],STDOUT_FILENO)==-1){// redirect stdout
            perror("redirecting stdout");
            return;
        }
        if(dup2(StderrPipe[PIPE_WRITE],STDERR_FILENO)==-1){// redirect stderr
            perror("redirecting stderr");
            return;
        }
        close(StdinPipe[PIPE_READ]);
        close(StdinPipe[PIPE_WRITE]);
        close(StdoutPipe[PIPE_READ]);
        close(StdoutPipe[PIPE_WRITE]);
        close(StderrPipe[PIPE_READ]);
        close(StderrPipe[PIPE_WRITE]);
        execl(Bot.name.c_str(),Bot.name.c_str(),(char*)NULL);//(char*)Null is really important
        //If you get past the previous line its an error
        perror("exec of the child process");
    }
    else if(nchild>0){//Parent process
        close(StdinPipe[PIPE_READ]);//Parent does not read from stdin of child
        close(StdoutPipe[PIPE_WRITE]);//Parent does not write to stdout of child
        close(StderrPipe[PIPE_WRITE]);//Parent does not write to stderr of child
        Bot.inPipe=StdinPipe[PIPE_WRITE];
        Bot.outPipe=StdoutPipe[PIPE_READ];
        Bot.errPipe=StderrPipe[PIPE_READ];
        Bot.pid=nchild;
    }
    else{//failed to create child
        close(StdinPipe[PIPE_READ]);
        close(StdinPipe[PIPE_WRITE]);
        close(StdoutPipe[PIPE_READ]);
        close(StdoutPipe[PIPE_WRITE]);
        perror("Failed to create child process");
    }
}

inline bool IsValidMove(const state &S,const AI &Bot,const string &M){
    return count(M.begin(),M.end(),'\n')==1;
}

string GetMove(const state &S,AI &Bot,const int turn){
    pollfd outpoll{Bot.outPipe,POLLIN};
    time_point<system_clock> Start_Time{system_clock::now()};
    string out;
    while(static_cast<duration<double>>(system_clock::now()-Start_Time).count()<(turn==1?FirstTurnTime:TimeLimit) && !IsValidMove(S,Bot,out)){
        double TimeLeft{(turn==1?FirstTurnTime:TimeLimit)-static_cast<duration<double>>(system_clock::now()-Start_Time).count()};
        if(poll(&outpoll,1,TimeLeft)){
            out+=EmptyPipe(Bot.outPipe);
        }
    }
    return out;
}

inline bool Has_Won(const array<AI,N> &Bot,const int idx)noexcept{
    if(!Bot[idx].alive()){
        return false;
    }
    for(int i=0;i<N;++i){
        if(i!=idx && Bot[i].alive()){
            return false;
        }
    }
    return true;
}

inline bool All_Dead(const array<AI,N> &Bot)noexcept{
    for(const AI &b:Bot){
        if(b.alive()){
            return false;
        }
    }
    return true;
}

action StringToAction(const state &S,const string &M_Str){
    action M;
    stringstream ss(M_Str);
    string type;
    ss >> type;
    if(type=="GOTO"){
        string destination;
        ss >> destination;
        if(StrToLocation.find(destination)==StrToLocation.end()){//Invalid destination
            cerr << "Invalid destination: " << M_Str << endl;
            throw(3);
        }
        M=action{GOTO,StrToLocation.at(destination)};
    }
    else if(type=="CONNECT"){
        string type;
        ss >> type;
        if(TypeToInt.find(type)!=TypeToInt.end()){
            M=action{CONNECT,TypeToInt.at(type)};
        }
        else{
            try{
                M=action{CONNECT,stoi(type)};
            }
            catch(...){//Invalid connect, neither a type nor an id
                cerr << "Invalid CONNECT: " << M_Str << endl;
                throw(3);
            }
        }
    }
    else{//Invalid move
        cerr << "Invalid move: " << M_Str << endl;
        throw(3);
    }
    return M;
}

inline bool ValidMoleculeIndex(const int idx)noexcept{
    return idx>=0 && idx<=4;
}

inline bool ValidRank(const int rank)noexcept{
    return rank>0 && rank<=3;
}

bool ReadyToProduce(const player &p,const sample &s){
    if(!s.diagnosed){
        return false;
    }
    for(int i=0;i<5;++i){//Look for missing molecules
        if(p.Mol[i]<s.Cost[i]-p.Exp[i]){
            return false;
        }
    }
    return true;
}

bool Completed_Project(const player &p,const project &proj){
    for(int m=0;m<5;++m){
        if(p.Exp[m]<proj.Target[m]){
            return false;
        }
    }
    return true;
}

void Simulate(state &S,const array<action,N> &M){
    const state S_before=S;
    for(int i=0;i<2;++i){
        player& p{S.P[i]};
        const action& mv{M[i]};
        if(p.eta==0){//Ignore actions of moving players
            if(mv.type==GOTO){
                p.eta=Distances[p.r][mv.id];
                p.r=intToLocation[mv.id];
            }
            else{//Connect
                if(p.r==SAMPLES && ValidRank(mv.id)){//Take undiagnosed sample
                    const int& rank{mv.id};
                    S.SamplePool[rank-1].push_back(S.SamplePool[rank-1].front());//Make a copy of the taken sample at the back of the list
                    sample s=S.SamplePool[rank-1].front();
                    S.SamplePool[rank-1].pop_front();
                    s.id=S.SampleCount++;
                    s.rank=rank;
                    s.diagnosed=false;
                    s.owner=i;
                    p.Samp.push_back(s);
                }
                else if(p.r==MOLECULES && ValidMoleculeIndex(mv.id)){//Take molecule
                    if(S_before.Avail[mv.id]>0 && accumulate(p.Mol.begin(),p.Mol.end(),0)<10){
                        --S.Avail[mv.id];
                        ++p.Mol[mv.id];
                    }
                }
                else if(p.r==LABORATORY && find_if(p.Samp.begin(),p.Samp.end(),[&](const sample &s){return s.id==mv.id;})!=p.Samp.end()){
                    const auto s=find_if(p.Samp.begin(),p.Samp.end(),[&](const sample &s){return s.id==mv.id;});
                    if(ReadyToProduce(p,*s)){
                        for(int m=0;m<5;++m){
                            const int spent{max(0,s->Cost[m]-p.Exp[m])};
                            p.Mol[m]-=spent;
                            S.Avail[m]+=spent;
                        }
                        p.score+=s->score;//Increase score
                        ++p.Exp[s->exp];//Gain expertise
                        p.Samp.erase(s);
                    }
                    else{
                        cerr << i << " tried to produce something he can't" << endl;
                    }
                }
                else if(p.r==DIAGNOSIS){
                    const auto player_s{find_if(p.Samp.begin(),p.Samp.end(),[&](const sample &s){return s.id==mv.id;})};
                    const auto diag_s{find_if(S.Samp.begin(),S.Samp.end(),[&](const sample &s){return s.id==mv.id;})};
                    if(player_s!=p.Samp.end()){
                        if(player_s->diagnosed){
                            S.Samp.push_back(*player_s);
                            p.Samp.erase(player_s);
                        }
                        else{
                            player_s->diagnosed=true;
                        }
                    }
                    else if(diag_s!=S.Samp.end()){
                        const action mv2{M[(i+1)%2]};
                        const player& p2{S.P[(i+1)%2]};
                        const bool other_hasnt_requested{mv.id!=mv2.id || mv.type!=mv2.type || p2.r!=DIAGNOSIS || p2.eta>0 || p2.Samp.size()==3};
                        if(p.Samp.size()<3 && (diag_s->owner==i || other_hasnt_requested) ){
                            //cerr << "Player " << i << " got sample " << mv.id << endl;
                            p.Samp.push_back(*diag_s);
                            S.Samp.erase(diag_s);
                        }
                    }
                }
            }
        }
    }
    for(int i=0;i<2;++i){//Decrease eta of both players
        player& p{S.P[i]};
        p.eta=max(0,p.eta-1);
    }
    for(auto it=S.Proj.begin();it!=S.Proj.end();){
        bool completed{false};
        for(int i=0;i<2;++i){
            if(Completed_Project(S.P[i],*it)){
                //cerr << "Player " << i << " completed a project" << endl;
                S.P[i].score+=50;
                completed=true;
            }
        }
        if(completed){
            it=S.Proj.erase(it);
        }
        else{
            ++it;
        }
    }
}

int Play_Game(const array<string,N> &Bot_Names,state &S){
    array<AI,N> Bot;
    for(int i=0;i<N;++i){
        Bot[i].id=i;
        Bot[i].name=Bot_Names[i];
        StartProcess(Bot[i]);
        stringstream ss;
        ss << S.Proj.size() << endl;
        for(int p=0;p<S.Proj.size();++p){//Project inputs
            for(int m=0;m<5;++m){
                ss << S.Proj[p].Target[m] << " ";
            }
            ss << endl;
        }
        Bot[i].Feed_Inputs(ss.str());
    }
    int turn{0};
    while(++turn>0 && !stop){
        array<action,N> M;
        for(int i=0;i<N;++i){
            if(Bot[i].alive()){
                stringstream ss;
                for(int j=0;j<2;++j){//Two player inputs
                    const player& p{S.P[(i+j)%2]};
                    ss << locationToStr[p.r] << " " << p.eta << " " << p.score;
                    for(int m=0;m<5;++m){
                        ss << " " << p.Mol[m];
                    }
                    for(int m=0;m<5;++m){
                        ss << " " << p.Exp[m];
                    }
                    ss << endl;
                }
                for(int m=0;m<5;++m){
                    ss << max(0,S.Avail[m]) << " ";
                }
                ss << endl;
                ss << S.Samp.size()+S.P[0].Samp.size()+S.P[1].Samp.size() << endl;
                for(const sample &s:S.Samp){
                    ss << s.id << " " << -1 << " " << s.rank << " " << typeToStr[s.exp] << " " << s.score;
                    for(int m=0;m<5;++m){
                        ss << " " << s.Cost[m];
                    }
                    ss << endl;
                }
                for(int j=0;j<2;++j){
                    const player& p{S.P[(i+j)%2]};
                    for(const sample &s:p.Samp){
                        ss << s.id << " " << j << " " << s.rank << " " << (s.diagnosed?typeToStr[s.exp]:"0") << " " << (s.diagnosed?s.score:-1);
                        for(int m=0;m<5;++m){
                            ss << " " << (s.diagnosed?s.Cost[m]:-1);
                        }
                        ss << endl;
                    }
                }
                try{
                    Bot[i].Feed_Inputs(ss.str());
                    string out=GetMove(S,Bot[i],turn);
                    //cerr << i << " " << out << endl;
                    M[i]=StringToAction(S,out);
                }
                catch(int ex){
                    if(ex==1){//Timeout
                        cerr << "Loss by Timeout of AI " << Bot[i].id << " name: " << Bot[i].name << endl;
                    }
                    else if(ex==3){
                        cerr << "Invalid move from AI " << Bot[i].id << " name: " << Bot[i].name << endl;
                    }
                    else if(ex==4){
                        cerr << "Error emptying pipe of AI " << Bot[i].name << endl;
                    }
                    else if(ex==5){
                        cerr << "AI " << Bot[i].name << " died before being able to give it inputs" << endl;
                    }
                    Bot[i].stop(turn);
                }
            }
        }
        Simulate(S,M);
        for(int i=0;i<N;++i){
            string err_str{EmptyPipe(Bot[i].errPipe)};
            if(Debug_AI){
                ofstream err_out("log.txt",ios::app);
                err_out << err_str << endl;
            }
        }
        for(int i=0;i<N;++i){
            if(Has_Won(Bot,i)){
                //cerr << i << " has won in " << turn << " turns" << endl;
                return i;
            }
        }
        if(All_Dead(Bot)){
            return -1;
        }
        if(turn==200){
            //cerr << S.P[0].score << " " << S.P[1].score << endl;
            return S.P[0].score==S.P[1].score?-1:S.P[0].score>S.P[1].score?0:1;
        }
    }
    return -2;
}

int Play_Round(array<string,N> Bot_Names){
    default_random_engine generator(system_clock::now().time_since_epoch().count());
    uniform_int_distribution<int> Swap_Distrib(0,1);
    const bool player_swap{Swap_Distrib(generator)==1};
    if(player_swap){
        swap(Bot_Names[0],Bot_Names[1]);
    }
    //Initial state generation
    state S;
    for(player &p:S.P){
        p=player{START,0,0,{0,0,0,0,0},{0,0,0,0,0},{}};
    }
    for(int &m:S.Avail){
        m=5;//5 molecules per type
    }
    array<vector<sample_template>,3> Pool_Vec=SampleList;
    for(int i=0;i<3;++i){
        random_shuffle(Pool_Vec[i].begin(),Pool_Vec[i].end());
        copy(Pool_Vec[i].begin(),Pool_Vec[i].end(),back_inserter(S.SamplePool[i]));
    }
    vector<project> Projects=ProjectList;
    random_shuffle(Projects.begin(),Projects.end());
    for(int i=0;i<3;++i){
        S.Proj.push_back(Projects[i]);
    }
    S.SampleCount=0;

    int winner{Play_Game(Bot_Names,S)};
    if(player_swap){
        return winner==-1?-1:winner==0?1:0;
    }
    else{
        return winner;
    }
}

void StopArena(const int signum){
    stop=true;
}

int main(int argc,char **argv){
    if(argc<3){
        cerr << "Program takes 2 inputs, the names of the AIs fighting each other" << endl;
        return 0;
    }
    int N_Threads{1};
    if(argc>=4){//Optional N_Threads parameter
        N_Threads=min(2*omp_get_num_procs(),max(1,atoi(argv[3])));
        cerr << "Running " << N_Threads << " arena threads" << endl;
    }
    array<string,N> Bot_Names;
    for(int i=0;i<2;++i){
        Bot_Names[i]=argv[i+1];
    }
    cout << "Testing AI " << Bot_Names[0];
    for(int i=1;i<N;++i){
        cerr << " vs " << Bot_Names[i];
    }
    cerr << endl;
    for(int i=0;i<N;++i){//Check that AI binaries are present
        ifstream Test{Bot_Names[i].c_str()};
        if(!Test){
            cerr << Bot_Names[i] << " couldn't be found" << endl;
            return 0;
        }
        Test.close();
    }
    signal(SIGTERM,StopArena);//Register SIGTERM signal handler so the arena can cleanup when you kill it
    signal(SIGPIPE,SIG_IGN);//Ignore SIGPIPE to avoid the arena crashing when an AI crashes
    int games{0},draws{0};
    array<double,2> points{0,0};
    #pragma omp parallel num_threads(N_Threads) shared(games,points,Bot_Names)
    while(!stop){
        int winner{Play_Round(Bot_Names)};
        if(winner==-1){//Draw
            #pragma omp atomic
            ++draws;
            #pragma omp atomic
            points[0]+=0.5;
            #pragma omp atomic
            points[1]+=0.5;
        }
        else{//Win
            ++points[winner];
        }
        #pragma omp atomic
        ++games;
        double p{static_cast<double>(points[0])/games};
        double sigma{sqrt(p*(1-p)/games)};
        double better{0.5+0.5*erf((p-0.5)/(sqrt(2)*sigma))};
        #pragma omp critical
        cout << "Wins:" << setprecision(4) << 100*p << "+-" << 100*sigma << "% Rounds:" << games << " Draws:" << draws << " " << better*100 << "% chance that " << Bot_Names[0] << " is better" << endl;
    }
}