/*
 *  /MathJax/jax/output/CommonHTML/autoload/multiline.js
 *
 *  Copyright (c) 2009-2018 The MathJax Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

MathJax.Hub.Register.StartupHook("CommonHTML Jax Ready",function(){var e="2.7.6";var b=MathJax.ElementJax.mml,a=MathJax.Hub.config,c=MathJax.OutputJax.CommonHTML;var g=b.mo().With({CHTML:c.BBOX.empty()});var f={newline:0,nobreak:1000000,goodbreak:[-200],badbreak:[+200],auto:[0],maxwidth:1.33,toobig:800,nestfactor:400,spacefactor:-100,spaceoffset:2,spacelimit:1,fence:500,close:500};var d={linebreakstyle:"after"};b.mbase.Augment({CHTMLlinebreakPenalty:f,CHTMLmultiline:function(k){var q=this;while(q.inferred||(q.parent&&q.parent.type==="mrow"&&q.parent.isEmbellished())){q=q.parent}var o=((q.type==="math"&&q.Get("display")==="block")||q.type==="mtd");q.isMultiline=true;var r=this.getValues("linebreak","linebreakstyle","lineleading","linebreakmultchar","indentalign","indentshift","indentalignfirst","indentshiftfirst","indentalignlast","indentshiftlast");if(r.linebreakstyle===b.LINEBREAKSTYLE.INFIXLINEBREAKSTYLE){r.linebreakstyle=this.Get("infixlinebreakstyle")}r.lineleading=this.CHTMLlength2em(r.lineleading,0.5);c.BBOX.empty(this.CHTML);var p=c.addElement(k,"mjx-stack");var h={BBOX:this.CHTML,n:0,Y:0,scale:(this.CHTML.scale||1),isTop:o,values:{},VALUES:r},n=this.CHTMLgetAlign(h,{}),j=this.CHTMLgetShift(h,{},n),i=[],l={index:[],penalty:f.nobreak,w:0,W:j,shift:j,scanW:j,nest:0},m=false;while(this.CHTMLbetterBreak(l,h,true)&&(l.scanW>=c.linebreakWidth||l.penalty===f.newline)){this.CHTMLaddLine(p,i,l.index,h,l.values,m);i=l.index.slice(0);m=true;n=this.CHTMLgetAlign(h,l.values);j=this.CHTMLgetShift(h,l.values,n);l.W=l.shift=l.scanW=j;l.penalty=f.nobreak}h.isLast=true;this.CHTMLaddLine(p,i,[],h,d,m);k.style.width=p.style.width=this.CHTML.pwidth="100%";this.CHTML.mwidth=c.Em(this.CHTML.w);this.CHTML.isMultiline=q.CHTML.isMultiline=true;p.style.verticalAlign=c.Em(h.d-this.CHTML.d);return k},CHTMLbetterBreak:function(l,h,s){if(this.isToken){return false}if(this.isEmbellished()){l.embellished=this;return this.CoreMO().CHTMLbetterBreak(l,h)}if(this.linebreakContainer){return false}var r=l.index.slice(0),p=l.index.shift(),o=this.data.length,n,t,k,q=(l.index.length>0),j=false;if(p==null){p=-1}if(!q){p++;l.W+=l.w;l.w=0}k=l.scanW=l.W;l.nest++;while(p<o&&(l.scanW<f.maxwidth*c.linebreakWidth||l.w===0)){if(this.data[p]){if(this.data[p].CHTMLbetterBreak(l,h)){j=true;r=[p].concat(l.index);n=l.W;t=l.w;if(l.penalty===f.newline){l.index=r;if(l.nest){l.nest--}return true}}k=(q?l.scanW:this.CHTMLaddWidth(p,l,k))}l.index=[];p++;q=false}if(s&&j){g.parent=this.parent;g.inherit=this.inherit;if(g.CHTMLbetterBreak(l,h)){j=false;r=l.index}}if(l.nest){l.nest--}l.index=r;if(j){l.W=n;l.w=t}return j},CHTMLaddWidth:function(h,k,j){if(this.data[h]){var l=this.data[h].CHTML;j+=(l.w+(l.L||0)+(l.R||0))*(l.scale||1);k.W=k.scanW=j;k.w=0}return j},CHTMLaddLine:function(s,j,o,h,t,q){var n=c.addElement(s,"mjx-block",{},[["mjx-box"]]),v=n.firstChild;var u=h.bbox=c.BBOX.empty();h.first=q;h.last=true;this.CHTMLmoveLine(j,o,v,h,t);u.clean();var r=this.CHTMLgetAlign(h,t),k=this.CHTMLgetShift(h,t,r,true);var m=0;if(h.n>0){var p=c.FONTDATA.baselineskip;var l=(h.values.lineleading==null?h.VALUES:h.values).lineleading*h.scale;var i=h.Y;h.Y-=Math.max(p,h.d+u.h+l);m=i-h.Y-h.d-u.h}if(k){v.style.margin="0 "+c.Em(-k)+" 0 "+c.Em(k)}if(r!==b.INDENTALIGN.LEFT){n.style.textAlign=r}if(m){n.style.paddingTop=c.Em(m)}h.BBOX.combine(u,k,h.Y);h.d=h.bbox.d;h.values=t;h.n++},CHTMLgetAlign:function(k,h){var l=h,i=k.values,j=k.VALUES,m;if(k.n===0){m=l.indentalignfirst||i.indentalignfirst||j.indentalignfirst}else{if(k.isLast){m=i.indentalignlast||j.indentalignlast}else{m=i.indentalign||j.indentalign}}if(m===b.INDENTALIGN.INDENTALIGN){m=i.indentalign||j.indentalign}if(m===b.INDENTALIGN.AUTO){m=(k.isTop?a.displayAlign:b.INDENTALIGN.LEFT)}return m},CHTMLgetShift:function(h,p,n,l){var o=p,m=h.values,i=h.VALUES,j;if(h.n===0){j=o.indentshiftfirst||m.indentshiftfirst||i.indentshiftfirst}else{if(h.isLast){j=m.indentshiftlast||i.indentshiftlast}else{j=m.indentshift||i.indentshift}}if(j===b.INDENTSHIFT.INDENTSHIFT){j=m.indentshift||i.indentshift}if(j==="auto"||j===""){j="0"}j=this.CHTMLlength2em(j,c.cwidth);if(h.isTop&&a.displayIndent!=="0"){var k=this.CHTMLlength2em(a.displayIndent,c.cwidth);j+=(n===b.INDENTALIGN.RIGHT?-k:k)}return(n===b.INDENTALIGN.RIGHT&&!l?-j:j)},CHTMLmoveLine:function(q,h,o,p,k){var m=q[0],l=h[0];if(m==null){m=-1}if(l==null){l=this.data.length-1}if(m===l&&q.length>1){this.data[m].CHTMLmoveSlice(q.slice(1),h.slice(1),o,p,k,"marginLeft")}else{var n=p.last;p.last=false;while(m<l){if(this.data[m]){if(q.length<=1){this.data[m].CHTMLmoveNode(o,p,k)}else{this.data[m].CHTMLmoveSlice(q.slice(1),[],o,p,k,"marginLeft")}}m++;p.first=false;q=[]}p.last=n;if(this.data[m]){if(h.length<=1){this.data[m].CHTMLmoveNode(o,p,k)}else{this.data[m].CHTMLmoveSlice([],h.slice(1),o,p,k,"marginRight")}}}},CHTMLmoveSlice:function(n,h,j,l,i,k){var m=this.CHTMLcreateSliceNode(j);this.CHTMLmoveLine(n,h,m,l,i);if(m.style[k]){m.style[k]=""}if(this.CHTML.L){if(k!=="marginLeft"){l.bbox.w+=this.CHTML.L}else{m.className=m.className.replace(/ MJXc-space\d/,"")}}if(this.CHTML.R&&k!=="marginRight"){l.bbox.w+=this.CHTML.R}if(h.length===0){j=this.CHTMLnodeElement();if(this.href){j=j.parentNode}j.parentNode.removeChild(j);j.nextMathJaxNode.id=j.id}return m},CHTMLcreateSliceNode:function(j){var i=this.CHTMLnodeElement(),l=0;if(this.href){i=i.parentNode}var h=i;while(h.nextMathJaxNode){h=h.nextMathJaxNode;l++}var k=i.cloneNode(false);h.nextMathJaxNode=k;k.nextMathJaxNode=null;k.id+="-MJX-Continue-"+l;return j.appendChild(k)},CHTMLmoveNode:function(h,k,i){if(!(k.first||k.last)||(k.first&&k.values.linebreakstyle===b.LINEBREAKSTYLE.BEFORE)||(k.last&&i.linebreakstyle===b.LINEBREAKSTYLE.AFTER)){var j=this.CHTMLnodeElement();if(this.href){j=j.parentNode}h.appendChild(j);if(this.CHTML.pwidth&&!h.style.width){h.style.width=this.CHTML.pwidth}if(k.last){j.style.marginRight=""}if(k.first||k.nextIsFirst){j.style.marginLeft="";this.CHTML.L=0;j.className=j.className.replace(/ MJXc-space\d/,"")}if(k.first&&this.CHTML.w===0){k.nextIsFirst=true}else{delete k.nextIsFirst}k.bbox.combine(this.CHTML,k.bbox.w,0)}}});b.mfenced.Augment({CHTMLbetterBreak:function(n,h){var v=n.index.slice(0),t=n.index.shift(),q=this.data.length,p,x,o,u=(n.index.length>0),l=false;if(t==null){t=-1}if(!u){t++;n.W+=n.w;n.w=0}o=n.scanW=n.W;n.nest++;if(!this.dataI){this.dataI=[];if(this.data.open){this.dataI.push("open")}if(q){this.dataI.push(0)}for(var s=1;s<q;s++){if(this.data["sep"+s]){this.dataI.push("sep"+s)}this.dataI.push(s)}if(this.data.close){this.dataI.push("close")}}q=this.dataI.length;while(t<q&&(n.scanW<f.maxwidth*c.linebreakWidth||n.w===0)){var r=this.dataI[t];if(this.data[r]){if(this.data[r].CHTMLbetterBreak(n,h)){l=true;v=[t].concat(n.index);p=n.W;x=n.w;if(n.penalty===f.newline){n.index=v;if(n.nest){n.nest--}return true}}o=(u?n.scanW:this.CHTMLaddWidth(t,n,o))}n.index=[];t++;u=false}if(n.nest){n.nest--}n.index=v;if(l){n.W=p;n.w=x}return l},CHTMLmoveLine:function(l,o,m,h,s){var q=l[0],p=o[0];if(q==null){q=-1}if(p==null){p=this.dataI.length-1}if(q===p&&l.length>1){this.data[this.dataI[q]].CHTMLmoveSlice(l.slice(1),o.slice(1),m,h,s,"marginLeft")}else{var r=h.last;h.last=false;var n=this.dataI[q];while(q<p){if(this.data[n]){if(l.length<=1){this.data[n].CHTMLmoveNode(m,h,s)}else{this.data[n].CHTMLmoveSlice(l.slice(1),[],m,h,s,"marginLeft")}}q++;n=this.dataI[q];h.first=false;l=[]}h.last=r;if(this.data[n]){if(o.length<=1){this.data[n].CHTMLmoveNode(m,h,s)}else{this.data[n].CHTMLmoveSlice([],o.slice(1),m,h,s,"marginRight")}}}}});b.msubsup.Augment({CHTMLbetterBreak:function(k,h){if(!this.data[this.base]){return false}var p=k.index.slice(0),n=k.index.shift(),m,q,l,o=(k.index.length>0),j=false;if(!o){k.W+=k.w;k.w=0}l=k.scanW=k.W;if(n==null){this.CHTML.baseW=this.data[this.base].CHTML.w;this.CHTML.dw=this.CHTML.w-this.CHTML.baseW}if(this.data[this.base].CHTMLbetterBreak(k,h)){j=true;p=[this.base].concat(k.index);m=k.W;q=k.w;if(k.penalty===f.newline){j=o=true}}if(!o){this.CHTMLaddWidth(this.base,k,l)}k.scanW+=this.CHTML.dw;k.W=k.scanW;k.index=[];if(j){k.W=m;k.w=q;k.index=p}return j},CHTMLmoveLine:function(j,n,m,i,r){if(this.data[this.base]){var k=c.addElement(m,"mjx-base");if(j.length>1){this.data[this.base].CHTMLmoveSlice(j.slice(1),n.slice(1),k,i,r,"marginLeft")}else{if(n.length<=1){this.data[this.base].CHTMLmoveNode(k,i,r)}else{this.data[this.base].CHTMLmoveSlice([],n.slice(1),k,i,r,"marginRight")}}}if(n.length===0){var l=this.CHTMLnodeElement(),p=c.getNode(l,"mjx-stack"),o=c.getNode(l,"mjx-sup"),h=c.getNode(l,"mjx-sub");if(p){m.appendChild(p)}else{if(o){m.appendChild(o)}else{if(h){m.appendChild(h)}}}var q=i.bbox.w,s;if(o){s=this.data[this.sup].CHTML;i.bbox.combine(s,q,s.Y)}if(h){s=this.data[this.sub].CHTML;i.bbox.combine(s,q,s.Y)}}}});b.mmultiscripts.Augment({CHTMLbetterBreak:function(l,i){if(!this.data[this.base]){return false}var p=l.index.slice(0);l.index.shift();var n,q,m,o=(l.index.length>0),k=false;if(!o){l.W+=l.w;l.w=0}l.scanW=l.W;var r=this.CHTML,j=this.data[this.base].CHTML;var h=r.w-j.w-(r.X||0);l.scanW+=r.X||0;m=l.scanW;if(this.data[this.base].CHTMLbetterBreak(l,i)){k=true;p=[this.base].concat(l.index);n=l.W;q=l.w;if(l.penalty===f.newline){k=o=true}}if(!o){this.CHTMLaddWidth(this.base,l,m)}l.scanW+=h;l.W=l.scanW;l.index=[];if(k){l.W=n;l.w=q;l.index=p}return k},CHTMLmoveLine:function(m,p,o,j,v){var n,i=this.CHTMLbbox,u;if(m.length<1){n=this.CHTMLnodeElement();var r=c.getNode(n,"mjx-prestack"),s=c.getNode(n,"mjx-presup"),l=c.getNode(n,"mjx-presub");if(r){o.appendChild(r)}else{if(s){o.appendChild(s)}else{if(l){o.appendChild(l)}}}u=j.bbox.w;if(s){j.bbox.combine(i.presup,u+i.presup.X,i.presup.Y)}if(l){j.bbox.combine(i.presub,u+i.presub.X,i.presub.Y)}}if(this.data[this.base]){var k=c.addElement(o,"mjx-base");if(m.length>1){this.data[this.base].CHTMLmoveSlice(m.slice(1),p.slice(1),k,j,v,"marginLeft")}else{if(p.length<=1){this.data[this.base].CHTMLmoveNode(k,j,v)}else{this.data[this.base].CHTMLmoveSlice([],p.slice(1),k,j,v,"marginRight")}}}if(p.length===0){n=this.CHTMLnodeElement();var t=c.getNode(n,"mjx-stack"),q=c.getNode(n,"mjx-sup"),h=c.getNode(n,"mjx-sub");if(t){o.appendChild(t)}else{if(q){o.appendChild(q)}else{if(h){o.appendChild(h)}}}u=j.bbox.w;if(q){j.bbox.combine(i.sup,u,i.sup.Y)}if(h){j.bbox.combine(i.sub,u,i.sub.Y)}}}});b.mo.Augment({CHTMLbetterBreak:function(j,h){if(j.values&&j.values.id===this.CHTMLnodeID){return false}var p=this.getValues("linebreak","linebreakstyle","lineleading","linebreakmultchar","indentalign","indentshift","indentalignfirst","indentshiftfirst","indentalignlast","indentshiftlast","texClass","fence");if(p.linebreakstyle===b.LINEBREAKSTYLE.INFIXLINEBREAKSTYLE){p.linebreakstyle=this.Get("infixlinebreakstyle")}if(p.texClass===b.TEXCLASS.OPEN){j.nest++}if(p.texClass===b.TEXCLASS.CLOSE&&j.nest){j.nest--}var k=j.scanW;delete j.embellished;var o=this.CHTML.w+(this.CHTML.L||0)+(this.CHTML.R||0);if(p.linebreakstyle===b.LINEBREAKSTYLE.AFTER){k+=o;o=0}if(k-j.shift===0&&p.linebreak!==b.LINEBREAK.NEWLINE){return false}var l=c.linebreakWidth-k;if(h.n===0&&(p.indentshiftfirst!==h.VALUES.indentshiftfirst||p.indentalignfirst!==h.VALUES.indentalignfirst)){var m=this.CHTMLgetAlign(h,p),i=this.CHTMLgetShift(h,p,m);l+=(j.shift-i)}var n=Math.floor(l/c.linebreakWidth*1000);if(n<0){n=f.toobig-3*n}if(p.fence){n+=f.fence}if((p.linebreakstyle===b.LINEBREAKSTYLE.AFTER&&p.texClass===b.TEXCLASS.OPEN)||p.texClass===b.TEXCLASS.CLOSE){n+=f.close}n+=j.nest*f.nestfactor;var q=f[p.linebreak||b.LINEBREAK.AUTO]||0;if(!MathJax.Object.isArray(q)){if(q||l>=0){n=q*j.nest}}else{n=Math.max(1,n+q[0]*j.nest)}if(n>=j.penalty){return false}j.penalty=n;j.values=p;j.W=k;j.w=o;p.lineleading=this.CHTMLlength2em(p.lineleading,h.VALUES.lineleading);p.id=this.CHTMLnodeID;return true}});b.mspace.Augment({CHTMLbetterBreak:function(i,h){if(i.values&&i.values.id===this.CHTMLnodeID){return false}var o=this.getValues("linebreak");var l=o.linebreak;if(!l||this.hasDimAttr()){l=b.LINEBREAK.AUTO}var j=i.scanW,n=this.CHTML.w+(this.CHTML.L||0)+(this.CHTML.R||0);if(j-i.shift===0){return false}var k=c.linebreakWidth-j;var m=Math.floor(k/c.linebreakWidth*1000);if(m<0){m=f.toobig-3*m}m+=i.nest*f.nestfactor;var p=f[l]||0;if(l===b.LINEBREAK.AUTO&&n>=f.spacelimit&&!this.mathbackground&&!this.background){p=[(n+f.spaceoffset)*f.spacefactor]}if(!MathJax.Object.isArray(p)){if(p||k>=0){m=p*i.nest}}else{m=Math.max(1,m+p[0]*i.nest)}if(m>=i.penalty){return false}i.penalty=m;i.values=o;i.W=j;i.w=n;o.lineleading=h.VALUES.lineleading;o.linebreakstyle="before";o.id=this.CHTMLnodeID;return true}});MathJax.Hub.Register.StartupHook("TeX mathchoice Ready",function(){b.TeXmathchoice.Augment({CHTMLbetterBreak:function(i,h){return this.Core().CHTMLbetterBreak(i,h)},CHTMLmoveLine:function(l,h,j,k,i){return this.Core().CHTMLmoveSlice(l,h,j,k,i)}})});b.maction.Augment({CHTMLbetterBreak:function(i,h){return this.Core().CHTMLbetterBreak(i,h)},CHTMLmoveLine:function(l,h,j,k,i){return this.Core().CHTMLmoveSlice(l,h,j,k,i)}});b.semantics.Augment({CHTMLbetterBreak:function(i,h){return(this.data[0]?this.data[0].CHTMLbetterBreak(i,h):false)},CHTMLmoveLine:function(l,h,j,k,i){return(this.data[0]?this.data[0].CHTMLmoveSlice(l,h,j,k,i):null)}});MathJax.Hub.Startup.signal.Post("CommonHTML multiline Ready");MathJax.Ajax.loadComplete(c.autoloadDir+"/multiline.js")});
