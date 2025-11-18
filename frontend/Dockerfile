# Stage 1: Build the frontend
FROM node:20-alpine AS builder

WORKDIR /app

COPY package.json ./
COPY package-lock.json ./

RUN npm install --legacy-peer-deps
COPY . .

RUN npm run build

FROM caddy:alpine

COPY --from=builder /app/dist /usr/share/caddy


COPY Caddyfile /etc/caddy/Caddyfile

EXPOSE 80
