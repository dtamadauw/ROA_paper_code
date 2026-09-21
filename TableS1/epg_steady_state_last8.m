%TITLE: MATLAB code for "Motion-robust gradient-echo T2 mapping in the liver using readout alternation"
%Matlab & Python scripts to generate table S1 in the paper
%Author: Daiki Tamada
%Affiliation: Department of Radiology, University of Wisconsin-Madison
%Date: 9/21/2026
%Email: dtamada@wisc.edu

%By downloading, installing, or otherwise accessing or using the Software , you ("Recipient") agree to receive and use the above-
%identified SOFTWARE subject to the following terms, obligations and restrictions. If you do not agree to all of the following terms, 
%obligations and restrictions you are not permitted to download, install,
%execute, access, or use the SOFTWARE:

%1.	Originators of the SOFTWARE.  Provider is willing to license its rights in the SOFTWARE ("Provider's Rights") to academic researchers to use free of charge solely for academic, non-commercial research purposes subject to the terms and conditions outlined herein. The SOFTWARE was created at the University of Wisconsin ("UW") by Dakai Tamada. Please note Provider's Rights may include, but are not limited to, certain patents or patent applications owned by the Wisconsin Alumni Research Foundation ("WARF"). 
%2.	Limited License.  Provider hereby grants to Recipient a non-commercial, non-transferable, royalty-free, non-exclusive license, without the right to sublicense, under Provider's Rights to  download, install, access, execute and use the SOFTWARE solely for academic, non-commercial research purposes. SOFTWARE may not be used, directly or indirectly, to perform services for a fee or for the production or manufacture of products for sale to third parties. The foregoing license does not include any license to third party intellectual property that may be contained in the SOFTWARE; obtaining a license to such rights is Recipient's responsibility. 
%3.	Restrictions on SOFTWARE use and distribution.  Recipient shall not take, authorize, or permit any of the following actions with the SOFTWARE: (1) modify, translate or otherwise create any derivative works; or (2) publicly display (e.g., Internet) or publicly perform (e.g., present at a press conference); or (3) sell, lease, rent or lend; or (4) use it for any commercial purposes whatsoever. Recipient must fully reproduce and not obscure, alter or remove any of the Provider's proprietary notices that appear on the SOFTWARE, including copyright notices or additional license terms included with any the third party software contained in the SOFTWARE. Recipient may not provide any third party with access to the SOFTWARE or use the SOFTWARE on a timeshare or service bureau basis. Recipient represents that it is compliance with all applicable export control provisions and is not prohibited from receiving the SOFTWARE. 
%4.	Reservation of rights.  Provider retains all rights and title in the SOFTWARE, including without limitation all intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future exist in the SOFTWARE, regardless of form or medium. Provider retains ownership and all of Its rights in the SOFTWARE, including all of its intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future cover the SOFTWARE or any uses of the SOFTWARE, regardless of form or medium; title remains with Provider and the SOFTWARE is merely being loaned to Recipient for the specific purposes and under the specific restrictions stated herein. Nothing in this Agreement grants Recipient any additional rights to the SOFTWARE, any right to obtain any updates or new releases of the SOFTWARE, any commercial license for the SOFTWARE, or any other intellectual property owned or licensed by Provider. Provider has no obligation to provide any support, updates, or bug fixes.
%5.	Disclaimer of Warranty. PROVIDER IS PROVIDING THE SOFTWARE TO RECIPIENT ON AN "AS IS" BASIS. PROVIDER MAKES NO REPRESENTATIONS OR WARRANTIES CONCERNING THE SOFTWARE OR ANY OUTCOME THAT MAY BE OBTAINED BY USING THE SOFTWARE, AND EXPRESSLY DISCLAIMS ALL SUCH WARRANTIES, INCLUDING WITHOUT LIMITATION ANY EXPRESS OR IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT OF INTELLECTUAL PROPERTY RIGHTS. PROVIDER MAKES NO REMEDY THAT THE SOFTWARE WILL OPERATE ERROR FREE OR UNINTERRUPTED.
%6.	Limitation of Liability; Indemnity.  TO THE FULLEST EXTENT PERMITTED BY LAW, IN NO EVENT SHALL PROVIDER BE LIABLE TO RECIPIENT FOR ANY LOST PROFITS OR ANY DIRECT, INDIRECT, EXEMPLARY, PUNITIVE, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING FROM THE SOFTWARE OR ITS USE. FURTHERMORE, IN NO EVENT WILL PROVIDER'S LIABILITY TO RECIPIENT EXCEED $100. PROVIDER HAS NO LIABILITY FOR ANY DECISION, ACT OR OMISSION MADE BY RECIPIENT AS A RESULT OF USE OF THE SOFTWARE. To the extent permitted by applicable law, Recipient agrees to indemnify, defend and hold harmless Provider, UW, and the SOFTWARE authors against all claims and expenses, including legal expenses and reasonable attorneys fees, arising from Recipient's use of the SOFTWARE.
%7.	No use of names/trademarks.  Recipient shall not use Provider's name, or the name of any author of the SOFTWARE or that of UW, in any manner without the prior written approval of the entity or person whose name is being used.
%8.	Termination.  Without prejudice to any other rights, Provider may terminate this Agreement if Recipient fails to comply with the terms of this Agreement for any reason. Upon termination for any reason, Recipient must immediately destroy all copies of the SOFTWARE in Recipient's possession, custody, or control.	

function result = epg_steady_state_last8(p)
% Steady-state-initialized pathway EPG: 100 TR combined + final 8 TR traced.
if nargin<1, p=struct(); end
if ~isfield(p,'useROA'), error('p.useROA is required'); end
p=d(p,'NTR_total',100); p=d(p,'NTR_trace',8); p=d(p,'TR_s',15e-3);
p=d(p,'T1_s',1.2); p=d(p,'T2_s',0.1); p=d(p,'flip_deg',15);
p=d(p,'dkRO',1); p=d(p,'dkSL',1); p=d(p,'GRO_Tpm',1); p=d(p,'GSL_Tpm',0.5);
p=d(p,'ampPruneCombined',1e-14); p=d(p,'ampPruneTrace',1e-14);
p=d(p,'maxTracePaths',2e6); p=d(p,'kTol',1e-9);
if p.NTR_trace>=p.NTR_total, error('NTR_trace must be smaller than NTR_total'); end
E1=exp(-p.TR_s/p.T1_s); E2=exp(-p.TR_s/p.T2_s);
S=table(0,0,complex(0),complex(0),complex(1),'VariableNames',{'kRO','kSL','Fp','Fm','Z'});
sig=complex(zeros(p.NTR_total,1)); startTR=p.NTR_total-p.NTR_trace+1; S0=[];
for n=1:p.NTR_total
    if n==startTR, S0=S; end
    S=rf_combined(S,p.flip_deg);
    [dk,~]=grad_for_tr(p,n);
    S=evolve_combined(S,dk,E1,E2,p.ampPruneCombined);
    idx=abs(S.kRO)<=p.kTol & abs(S.kSL)<=p.kTol; sig(n)=sum(S.Fp(idx));
end
paths=states_to_paths(S0,p.ampPruneTrace);
for n=startTR:p.NTR_total
    paths=rf_branch(paths,p.flip_deg,n,p.ampPruneTrace);
    [dk,G]=grad_for_tr(p,n); t0=(n-1)*p.TR_s;
    paths=evolve_paths(paths,dk,G,t0,p.TR_s,E1,E2,n,p.ampPruneTrace);
    if numel(paths)>p.maxTracePaths
        [~,ii]=maxk(abs([paths.amp]),p.maxTracePaths); paths=paths(ii);
        warning('Trace paths capped at %d',p.maxTracePaths);
    end
end
keep=false(numel(paths),1);
for i=1:numel(paths), keep(i)=paths(i).state==1 && norm(paths(i).k)<=p.kTol; end
T=to_table(paths(keep));
result.parameters=p; result.paths=T; result.summary=summarize(T);
result.signalHistory=table((1:p.NTR_total).',sig,abs(sig),angle(sig),'VariableNames',{'TR','Signal','Magnitude','Phase_rad'});
result.decompositionCheck=table(sig(end),sum(T.Amplitude),abs(sig(end)-sum(T.Amplitude)),abs(sum(T.Amplitude))/max(abs(sig(end)),eps),...
    'VariableNames',{'CombinedSignal','DecomposedSignal','AbsoluteDifference','MagnitudeRatio'});
end

function S=rf_combined(S,aDeg)
R=rfmat(aDeg); for i=1:height(S), y=R*[S.Fp(i);S.Fm(i);S.Z(i)]; S.Fp(i)=y(1); S.Fm(i)=y(2); S.Z(i)=y(3); end
end
function S2=evolve_combined(S,dk,E1,E2,pr)
rows=repmat(struct('kRO',0,'kSL',0,'Fp',0,'Fm',0,'Z',0),0,1);
for i=1:height(S)
 if abs(S.Fp(i))>=pr, r=zrow(); r.kRO=S.kRO(i)+dk(1); r.kSL=S.kSL(i)+dk(2); r.Fp=S.Fp(i)*E2; rows(end+1)=r; end
 if abs(S.Fm(i))>=pr, r=zrow(); r.kRO=S.kRO(i)-dk(1); r.kSL=S.kSL(i)-dk(2); r.Fm=S.Fm(i)*E2; rows(end+1)=r; end
 if abs(S.Z(i))>=pr, r=zrow(); r.kRO=S.kRO(i); r.kSL=S.kSL(i); r.Z=S.Z(i)*E1; rows(end+1)=r; end
end
r=zrow(); r.Z=1-E1; rows(end+1)=r;
k=[[rows.kRO].',[rows.kSL].']; [ku,~,ic]=unique(k,'rows');
Fp=accumarray(ic,[rows.Fp].',[],@sum,0); Fm=accumarray(ic,[rows.Fm].',[],@sum,0); Z=accumarray(ic,[rows.Z].',[],@sum,0);
kk=abs(Fp)+abs(Fm)+abs(Z)>=pr; S2=table(ku(kk,1),ku(kk,2),Fp(kk),Fm(kk),Z(kk),'VariableNames',{'kRO','kSL','Fp','Fm','Z'});
end
function paths=states_to_paths(S,pr)
paths=repmat(newp(),0,1); states=[1,-1,0];
for i=1:height(S)
 vals=[S.Fp(i),S.Fm(i),S.Z(i)];
 for j=1:3
  if abs(vals(j))<pr, continue; end
  q=newp(); q.state=states(j); q.amp=vals(j); q.k=[S.kRO(i),S.kSL(i)];
  q.history=sprintf('SS(kRO=%g,kSL=%g):%s',S.kRO(i),S.kSL(i),sname(states(j)));
  q.startedTransverse=states(j)~=0; paths(end+1)=q;
 end
end
end
function out=rf_branch(paths,aDeg,n,pr)
R=rfmat(aDeg); out=repmat(newp(),0,1);
for i=1:numel(paths)
 col=statecol(paths(i).state);
 for row=1:3
  amp=paths(i).amp*R(row,col); if abs(amp)<pr, continue; end
  q=paths(i); q.state=rowstate(row); q.amp=amp;
  if q.state~=0, q.startedTransverse=true; elseif q.startedTransverse, q.returnedToZ=true; end
  q.history=q.history+">RF"+string(n)+":"+sname(q.state); out(end+1)=q;
 end
end
end
function out=evolve_paths(paths,dk,G,t0,dt,E1,E2,n,pr)
A1=G*(t0*dt+0.5*dt^2); out=repmat(newp(),0,1);
for i=1:numel(paths)
 q=paths(i);
 if q.state==0, q.amp=q.amp*E1; else, q.amp=q.amp*E2; q.k=q.k+q.state*dk; q.M1=q.M1+q.state*A1; end
 if abs(q.amp)>=pr, q.history=q.history+">TR"+string(n); out(end+1)=q; end
end
q=newp(); q.state=0; q.amp=1-E1; q.history="Recovery@TR"+string(n); out(end+1)=q;
end
function [dk,G]=grad_for_tr(p,n)
if p.useROA, s=(-1)^(n-1); else, s=1; end
dk=[s*p.dkRO,p.dkSL]; G=[s*p.GRO_Tpm,p.GSL_Tpm];
end
function T=to_table(paths)
n=numel(paths); Pathway=strings(n,1); EchoType=strings(n,1); Amplitude=complex(zeros(n,1)); Magnitude=zeros(n,1); Power=zeros(n,1); Phase_rad=zeros(n,1); M1_RO=zeros(n,1); M1_SL=zeros(n,1); M1_RSS=zeros(n,1);
for i=1:n
 q=paths(i); Pathway(i)=q.history; if q.returnedToZ, EchoType(i)="STE-like"; else, EchoType(i)="SE-like"; end
 Amplitude(i)=q.amp; Magnitude(i)=abs(q.amp); Power(i)=abs(q.amp)^2; Phase_rad(i)=angle(q.amp); M1_RO(i)=q.M1(1); M1_SL(i)=q.M1(2); M1_RSS(i)=norm(q.M1);
end
T=table(Pathway,EchoType,Amplitude,Magnitude,Power,Phase_rad,M1_RO,M1_SL,M1_RSS);
if isempty(T), T.RelativeMagnitude=zeros(0,1); T.RelativePower=zeros(0,1); return; end
T=sortrows(T,'Magnitude','descend'); T.RelativeMagnitude=T.Magnitude/sum(T.Magnitude); T.RelativePower=T.Power/sum(T.Power);
end
function S=summarize(T)
S=struct(); S.nPathways=height(T); if isempty(T), return; end
S.totalSignal=sum(T.Amplitude); w=T.Power; M=[T.M1_RO,T.M1_SL]; S.M1RMS_powerWeighted=sqrt(sum(w.*M.^2,1)/sum(w)); S.M1RSS_RMS_powerWeighted=sqrt(sum(w.*T.M1_RSS.^2)/sum(w));
is=T.EchoType=="STE-like"; S.fractionMagnitude_STE=sum(T.Magnitude(is))/sum(T.Magnitude); S.fractionPower_STE=sum(T.Power(is))/sum(T.Power);
end
function R=rfmat(aDeg)
a=deg2rad(aDeg); c2=cos(a/2)^2; s2=sin(a/2)^2; sa=sin(a); ca=cos(a); R=[c2,s2,-1i*sa;s2,c2,1i*sa;-0.5i*sa,0.5i*sa,ca];
end
function q=newp(), q=struct('state',0,'amp',0,'k',[0,0],'M1',[0,0],'history',"",'startedTransverse',false,'returnedToZ',false); end
function r=zrow(), r=struct('kRO',0,'kSL',0,'Fp',0,'Fm',0,'Z',0); end
function c=statecol(s), if s==1,c=1;elseif s==-1,c=2;else,c=3;end,end
function s=rowstate(r), a=[1,-1,0]; s=a(r); end
function s=sname(x), if x==1,s="F+";elseif x==-1,s="F-";else,s="Z";end,end
function p=d(p,n,v), if ~isfield(p,n)||isempty(p.(n)),p.(n)=v;end,end
